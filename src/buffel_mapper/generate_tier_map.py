#!/usr/bin/env python3
import os
# Optimizations for GDAL to stream COGs efficiently from cloud storage
os.environ["GDAL_HTTP_TIMEOUT"] = "60"
os.environ["GDAL_HTTP_CONNECTTIMEOUT"] = "30"
os.environ["GDAL_HTTP_MAX_RETRY"] = "5"
os.environ["GDAL_DISABLE_READDIR_ON_OPEN"] = "EMPTY_DIR"
os.environ["CPL_VSIL_CURL_ALLOWED_EXTENSIONS"] = ".tif"
os.environ["GDAL_HTTP_MERGE_CONSECUTIVE_RANGES"] = "YES"
os.environ["GDAL_CACHEMAX"] = "512"
os.environ["VSI_CACHE"] = "TRUE"
os.environ["VSI_CACHE_SIZE"] = "536870912"

import argparse
import numpy as np
import pandas as pd
import geopandas as gpd
import json
import pystac_client
import planetary_computer
import xarray as xr
import rasterio
import rioxarray
import shapely.geometry

# Force non-interactive matplotlib backend before importing pyplot to prevent interpreter shutdown errors
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Patch requests to have a default timeout of 15 seconds to prevent hanging
import requests
from urllib3.util import Retry
from requests.adapters import HTTPAdapter
session = requests.Session()
retries = Retry(total=3, backoff_factor=1, status_forcelist=[500, 502, 503, 504])
session.mount('http://', HTTPAdapter(max_retries=retries))
session.mount('https://', HTTPAdapter(max_retries=retries))

def patched_get(*args, **kwargs):
    if 'timeout' not in kwargs:
        kwargs['timeout'] = 15
    return session.get(*args, **kwargs)

def patched_post(*args, **kwargs):
    if 'timeout' not in kwargs:
        kwargs['timeout'] = 15
    return session.post(*args, **kwargs)

requests.get = patched_get
requests.post = patched_post

def rgb_to_float(rgb_array, cmap_name, vmin, vmax):
    h, w = rgb_array.shape[1], rgb_array.shape[2]
    rgb_3d = np.transpose(rgb_array, (1, 2, 0))
    flat_rgb = rgb_3d.reshape(-1, 3)
    is_bg = np.all(flat_rgb < 0.01, axis=1) | np.any(np.isnan(flat_rgb), axis=1)
    
    x = np.linspace(0, 1, 1000)
    cmap = plt.get_cmap(cmap_name)
    lut = cmap(x)[..., :3]
    
    raw_vals = np.full(flat_rgb.shape[0], np.nan)
    non_bg_rgb = flat_rgb[~is_bg]
    
    if len(non_bg_rgb) > 0:
        m_best_dist = np.full(non_bg_rgb.shape[0], np.inf)
        m_best_idx = np.zeros(non_bg_rgb.shape[0], dtype=int)
        for idx in range(1000):
            color = lut[idx]
            dist = np.sum((non_bg_rgb - color) ** 2, axis=1)
            closer = dist < m_best_dist
            m_best_dist[closer] = dist[closer]
            m_best_idx[closer] = idx
            
        norm_vals = m_best_idx / 999.0
        raw_vals[~is_bg] = norm_vals * (vmax - vmin) + vmin
        
    return raw_vals.reshape(h, w)

def load_raster_float(da_rio, cmap_name, vmin, vmax):
    """Read a single-band float32 raster directly, or decode a legacy 3-band RGB raster."""
    vals = da_rio.values
    if vals.shape[0] == 1:
        raw = vals[0].astype(float)
        nodata = da_rio.rio.nodata
        if nodata is not None and not np.isnan(nodata):
            raw[raw == nodata] = np.nan
        return raw
    return rgb_to_float(vals, cmap_name, vmin, vmax)

def dilate_mask_2d(mask_2d, radius=2):
    """Dilate a 2D boolean mask by `radius` pixels using pure NumPy."""
    if radius <= 0:
        return mask_2d.astype(bool)
    padded = np.pad(mask_2d.astype(bool), radius, mode='edge')
    h, w = mask_2d.shape
    out = np.zeros((h, w), dtype=bool)
    for dy in range(-radius, radius + 1):
        for dx in range(-radius, radius + 1):
            out |= padded[radius + dy : radius + dy + h, radius + dx : radius + dx + w]
    return out

def main():
    parser = argparse.ArgumentParser(description="Create 3-Tier Weed Management Spatial Layer")
    parser.add_argument("--project", type=str, default="Keewong_BSA_Shapefile_Resilient2", help="Project name")
    parser.add_argument("--buffer-pixels", type=int, default=2, help="Buffer pixels for paddock boundaries")
    
    args = parser.parse_args()
    
    pipeline_dir = os.path.dirname(os.path.abspath(__file__))
    if os.path.isdir(args.project):
        project_dir = os.path.abspath(args.project)
    elif os.path.isdir(os.path.join(os.getcwd(), "projects", args.project)):
        project_dir = os.path.join(os.getcwd(), "projects", args.project)
    else:
        project_dir = os.path.join(pipeline_dir, "projects", args.project)
    shapefile_path = os.path.join(project_dir, "BSA.shp")
    
    if not os.path.exists(shapefile_path):
        raise FileNotFoundError(f"Project shapefile not found at: {shapefile_path}")
        
    out_dir = os.path.join(project_dir, "tiers_layer")
    os.makedirs(out_dir, exist_ok=True)
    
    print("=" * 60)
    print(f"[*] CREATING 3-TIER MANAGEMENT LAYER FOR PROJECT: {args.project}")
    print("=" * 60)
    
    # 1. Load Paddock Geometry and Bounding Box
    gdf = gpd.read_file(shapefile_path)
    gdf_wgs84 = gdf.to_crs("EPSG:4326")
    # Dissolve to handle multi-row paddock geometries robustly
    geom_multi = gdf_wgs84.geometry.unary_union
    bbox_coords = list(geom_multi.bounds)
    
    # 2. Load Targeted Rasters (native float32 or legacy RGB fallback)
    max_spike_path = os.path.join(project_dir, "targeted", "buffel_max_growth_spike.tif")
    freq_path = os.path.join(project_dir, "targeted", "buffel_growth_frequency.tif")
    
    da_spike_rgb = rioxarray.open_rasterio(max_spike_path)
    da_freq_rgb = rioxarray.open_rasterio(freq_path)
    
    print("[*] Loading growth spike values...")
    spike_raw = load_raster_float(da_spike_rgb, 'YlGn', 0.0, 0.40)
    print("[*] Loading growth frequency values...")
    freq_raw = load_raster_float(da_freq_rgb, 'plasma', 0.0, 1.0)
    
    # 3. Connect to STAC API to compute 2025 mean and std MSAVI (background greenness/volatility)
    print("[*] Connecting to Earth Search STAC API...")
    catalog = pystac_client.Client.open(
        "https://earth-search.aws.element84.com/v1"
    )
    
    # Limit cloud cover to < 10% to speed up download and ensure clean statistics
    search = catalog.search(
        collections=["sentinel-2-l2a"],
        bbox=bbox_coords,
        datetime="2023-08-01/2025-11-15",
        query={"eo:cloud_cover": {"lt": 10}}
    )
    
    import time
    max_stac_retries = 3
    items = None
    for attempt in range(max_stac_retries):
        try:
            items = search.item_collection()
            break
        except Exception as e:
            if attempt == max_stac_retries - 1:
                raise e
            print(f"    [!] STAC query failed ({e}). Retrying in 5 seconds...", flush=True)
            time.sleep(5)
            
    # Filter items to only include dry season months (August, September, October, November)
    dry_season_items = []
    for it in items:
        dt_str = it.properties["datetime"]
        month = int(dt_str[5:7])
        if month in [8, 9, 10, 11]:
            dry_season_items.append(it)
            
    print(f"[+] Found {len(dry_season_items)} clean dry-season scenes across 2023-2025.")
    
    # Select cleanest 6 scenes across 2023-2025 for long-term Refugia baseline
    sorted_baseline_items = sorted(dry_season_items, key=lambda x: x.properties["eo:cloud_cover"])[:6]
    
    # Select earliest and latest scenes of the 2025 dry season for phenological curing rate
    items_2025 = [it for it in dry_season_items if "2025" in it.properties["datetime"]]
    items_2025_sorted = sorted(items_2025, key=lambda x: x.properties["datetime"])
    if len(items_2025_sorted) >= 2:
        earliest_2025 = items_2025_sorted[0]
        latest_2025 = items_2025_sorted[-1]
        curing_candidates = [earliest_2025, latest_2025]
    else:
        earliest_2025 = None
        latest_2025 = None
        curing_candidates = []
        
    # Combine into a list of unique items to download
    unique_items = []
    seen_ids = set()
    for item in sorted_baseline_items + curing_candidates:
        if item.id not in seen_ids:
            unique_items.append(item)
            seen_ids.add(item.id)
            
    msavi_list = []
    ndre_list = []
    ndwi_list = []
    crs_ref = da_spike_rgb.rio.crs
    
    # Create template DataArray matching the targeted raster grid shape
    da_spike_template = xr.DataArray(
        spike_raw,
        coords={'y': da_spike_rgb.y, 'x': da_spike_rgb.x},
        dims=['y', 'x']
    ).rio.write_crs(crs_ref)
    
    msavi_dict = {}
    ndre_dict = {}
    ndwi_dict = {}
    
    print(f"[*] Loading MSAVI, NDRE, and NDWI bands for {len(unique_items)} unique scenes...")
    
    for idx, item in enumerate(unique_items):
        try:
            date_str = item.properties["datetime"][:10]
            cloud_cover = item.properties["eo:cloud_cover"]
            print(f"    [*] Loading scene {idx+1}/{len(unique_items)} ({date_str}, cloud cover: {cloud_cover:.1f}%)...", flush=True)
            red_url = item.assets["red"].href
            nir_url = item.assets["nir"].href
            scl_url = item.assets["scl"].href
            b05_url = item.assets["rededge1"].href
            b11_url = item.assets["swir16"].href
            
            # Open Sentinel-2 bands (lazy metadata open)
            da_red = rioxarray.open_rasterio(red_url)
            da_nir = rioxarray.open_rasterio(nir_url)
            da_scl = rioxarray.open_rasterio(scl_url)
            da_b05 = rioxarray.open_rasterio(b05_url)
            da_b11 = rioxarray.open_rasterio(b11_url)
            
            # Reproject geometry once per scene to the scene's native CRS to avoid coordinate mismatches
            s2_crs = da_red.rio.crs
            gdf_s2 = gdf.to_crs(s2_crs)
            s2_bbox = gdf_s2.total_bounds
            geom_s2 = gdf_s2.geometry.unary_union.buffer(args.buffer_pixels * 10.0)
            
            # Clip bands to paddock geometry in the native scene CRS using fast bounding box pre-crop
            da_red_cropped = da_red.rio.clip_box(*s2_bbox).rio.clip([geom_s2], crs=s2_crs, all_touched=True)
            da_nir_cropped = da_nir.rio.clip_box(*s2_bbox).rio.clip([geom_s2], crs=s2_crs, all_touched=True)
            da_scl_cropped = da_scl.rio.clip_box(*s2_bbox).rio.clip([geom_s2], crs=s2_crs, all_touched=True)
            da_b05_cropped = da_b05.rio.clip_box(*s2_bbox).rio.clip([geom_s2], crs=s2_crs, all_touched=True)
            da_b11_cropped = da_b11.rio.clip_box(*s2_bbox).rio.clip([geom_s2], crs=s2_crs, all_touched=True)
            
            # Extract values, squeezing only the band dimension of size 1
            red_vals = da_red_cropped.values.astype(float).squeeze(axis=0) / 10000.0
            nir_vals = da_nir_cropped.values.astype(float).squeeze(axis=0) / 10000.0
            
            # Resample 20m SCL to 10m grid cleanly using rioxarray's reproject_match (nearest-neighbor)
            da_scl_resampled = da_scl_cropped.rio.reproject_match(
                da_red_cropped,
                resampling=rasterio.enums.Resampling.nearest
            )
            scl_vals = da_scl_resampled.values.squeeze(axis=0)
            
            # Resample 20m B05 and B11 to 10m grid using bilinear interpolation
            da_b05_resampled = da_b05_cropped.rio.reproject_match(
                da_red_cropped,
                resampling=rasterio.enums.Resampling.bilinear
            )
            da_b11_resampled = da_b11_cropped.rio.reproject_match(
                da_red_cropped,
                resampling=rasterio.enums.Resampling.bilinear
            )
            b05_vals = da_b05_resampled.values.astype(float).squeeze(axis=0) / 10000.0
            b11_vals = da_b11_resampled.values.astype(float).squeeze(axis=0) / 10000.0
            
            # Mask clouds and bad data (SCL 0=no_data, 1=saturated/defective, 3=cloud shadow, 8=med cloud, 9=high cloud, 10=cirrus)
            cloud_mask = dilate_mask_2d(np.isin(scl_vals, [0, 1, 3, 8, 9, 10]), radius=2)
            red_vals[cloud_mask] = np.nan
            nir_vals[cloud_mask] = np.nan
            b05_vals[cloud_mask] = np.nan
            b11_vals[cloud_mask] = np.nan
            
            # MSAVI formula
            sqrt_term = (2.0 * nir_vals + 1.0)**2 - 8.0 * (nir_vals - red_vals)
            sqrt_term = np.clip(sqrt_term, 0, None)  # Avoid sqrt of negative values from atmospheric correction offsets
            msavi = (2.0 * nir_vals + 1.0 - np.sqrt(sqrt_term)) / 2.0
            
            # NDRE formula (Red Edge chlorophyll sensitivity)
            ndre = (nir_vals - b05_vals) / np.clip(nir_vals + b05_vals, 1e-5, None)
            
            # NDWI formula (SWIR canopy moisture sensitivity)
            ndwi = (nir_vals - b11_vals) / np.clip(nir_vals + b11_vals, 1e-5, None)
            
            # Create 2D DataArrays in the native Sentinel-2 CRS
            da_msavi = xr.DataArray(
                msavi, 
                coords={'y': da_red_cropped.y, 'x': da_red_cropped.x}, 
                dims=['y', 'x']
            ).rio.write_crs(s2_crs)
            da_ndre = xr.DataArray(
                ndre,
                coords={'y': da_red_cropped.y, 'x': da_red_cropped.x},
                dims=['y', 'x']
            ).rio.write_crs(s2_crs)
            da_ndwi = xr.DataArray(
                ndwi,
                coords={'y': da_red_cropped.y, 'x': da_red_cropped.x},
                dims=['y', 'x']
            ).rio.write_crs(s2_crs)
            
            # Align immediately to target template grid to ensure coordinate alignment before concat
            da_msavi_aligned = da_msavi.rio.reproject_match(da_spike_template)
            da_ndre_aligned = da_ndre.rio.reproject_match(da_spike_template)
            da_ndwi_aligned = da_ndwi.rio.reproject_match(da_spike_template)
            
            msavi_dict[item.id] = da_msavi_aligned
            ndre_dict[item.id] = da_ndre_aligned
            ndwi_dict[item.id] = da_ndwi_aligned
            print(f"    [+] Loaded scene {idx+1}/{len(unique_items)}", flush=True)
        except Exception as e:
            import traceback
            traceback.print_exc()
            print(f"    [!] Error loading scene {idx+1}: {e}", flush=True)
            continue
            
    loaded_baseline_items = [it for it in sorted_baseline_items if it.id in msavi_dict]
    if len(loaded_baseline_items) < 2:
        raise RuntimeError("Too few successful baseline scenes loaded to build a reliable baseline.")
        
    print("[*] Computing mean, minimum, and standard deviation of MSAVI, NDRE, and NDWI dry-season baseline...")
    msavi_stack = xr.concat([msavi_dict[it.id] for it in loaded_baseline_items], dim='time')
    mean_msavi_aligned = msavi_stack.mean(dim='time', skipna=True).values.squeeze()
    min_msavi_aligned = msavi_stack.min(dim='time', skipna=True).values.squeeze()
    std_msavi_aligned = msavi_stack.std(dim='time', skipna=True).values.squeeze()
    
    ndre_stack = xr.concat([ndre_dict[it.id] for it in loaded_baseline_items], dim='time')
    mean_ndre_aligned = ndre_stack.mean(dim='time', skipna=True).values.squeeze()
    
    ndwi_stack = xr.concat([ndwi_dict[it.id] for it in loaded_baseline_items], dim='time')
    mean_ndwi_aligned = ndwi_stack.mean(dim='time', skipna=True).values.squeeze()
    
    # Dry-season senescence verification: confirms grass cures down during dry spells (unlike permanent evergreen canopy)
    print("[*] Computing dry-season senescence mask (min dry-season MSAVI < 0.22)...")
    curing_mask = np.isnan(min_msavi_aligned) | (min_msavi_aligned < 0.22)
        
    # 5. Apply the 3-Tier Decision Tree Logic
    print("[*] Applying Tier Delineation logic...")
    # Initialize tier layer (0: Background/Soil, 1: Tier 1, 2: Tier 2, 3: Tier 3)
    tiers = np.zeros_like(spike_raw, dtype=np.uint8)
    
    # Valid data mask (non-NaN in targeted growth spike)
    valid_mask = ~np.isnan(spike_raw)
    valid_msavi = ~np.isnan(mean_msavi_aligned) & ~np.isnan(std_msavi_aligned) & ~np.isnan(mean_ndre_aligned) & ~np.isnan(mean_ndwi_aligned)
    
    # Calibrate thresholds based on 90th-percentile MSAVI (1-band float32) vs legacy 3-band Max NDVI
    is_p90_msavi = (da_spike_rgb.shape[0] == 1)
    t3_spike_max = 0.085 if is_p90_msavi else 0.12
    t2_spike_min = 0.105 if is_p90_msavi else 0.20
    t2_freq_min  = 0.050 if is_p90_msavi else 0.30
    t1_spike_min = 0.080 if is_p90_msavi else 0.12
    t1_freq_min  = 0.025 if is_p90_msavi else 0.10

    # Tier 3: Native Refugia (Woody Evergreen Canopy Enhanced)
    # - Consistently green: mean MSAVI >= 0.20
    # - Active woody chlorophyll: mean NDRE >= 0.15 (Red Edge sensitivity)
    # - Retains leaf moisture in drought: mean NDWI >= -0.22 (SWIR sensitivity, calibrated for dryland soil backgrounds)
    # - Low volatility: standard deviation < 0.06
    # - Stable: 90th-percentile growth spike < t3_spike_max (does not respond to rain pulses like grass)
    tier3_mask = (
        valid_mask & valid_msavi & 
        (mean_msavi_aligned >= 0.20) & 
        (mean_ndre_aligned >= 0.15) & 
        (mean_ndwi_aligned >= -0.22) & 
        (std_msavi_aligned < 0.06) & 
        (spike_raw < t3_spike_max)
    )
    
    # Tier 2: Heavy Infestations (Core weed zones)
    tier2_mask = valid_mask & (spike_raw >= t2_spike_min) & (freq_raw >= t2_freq_min) & (~tier3_mask) & curing_mask
    
    # Tier 1: Light & New Infestations (Active eradication zones)
    tier1_mask = valid_mask & ((spike_raw >= t1_spike_min) & (freq_raw >= t1_freq_min)) & (~tier3_mask) & (~tier2_mask) & curing_mask
    
    # Assign values
    tiers[tier1_mask] = 1
    tiers[tier2_mask] = 2
    tiers[tier3_mask] = 3
    
    # Mask out background areas (NaN outside KML)
    tiers_float = tiers.astype(float)
    tiers_float[~valid_mask] = np.nan
    
    # 6. Save outputs
    transform_ref = da_spike_rgb.rio.transform()
    
    # Save raw single-band GeoTIFF
    out_tif_path = os.path.join(out_dir, "buffel_management_tiers.tif")
    meta = {
        'driver': 'GTiff',
        'dtype': 'float32',
        'nodata': -9999.0,
        'width': tiers.shape[1],
        'height': tiers.shape[0],
        'count': 1,
        'crs': crs_ref,
        'transform': transform_ref
    }
    
    tiers_output = tiers_float.copy()
    tiers_output[np.isnan(tiers_output)] = -9999.0
    with rasterio.open(out_tif_path, 'w', **meta) as dst:
        dst.write(tiers_output, 1)
    print(f"[+] Saved management tiers GeoTIFF to: {out_tif_path}")
    
    # Reproject tiers to EPSG:4326 for Leaflet visualization to avoid spatial distortion/shift
    da_tiers = xr.DataArray(
        tiers_float,
        coords={'y': da_spike_rgb.y, 'x': da_spike_rgb.x},
        dims=['y', 'x']
    ).rio.write_crs(crs_ref)
    
    da_tiers_gps = da_tiers.rio.reproject("EPSG:4326", resampling=rasterio.enums.Resampling.nearest)
    tiers_gps = da_tiers_gps.values.squeeze()
    
    # Colors: 0: Transparent/gray, 1: Green (#22c55e), 2: Red (#ef4444), 3: Cyan (#06b6d4)
    # We create a custom discrete colormap
    from matplotlib.colors import ListedColormap
    discrete_cmap = ListedColormap(['#4b5563', '#22c55e', '#ef4444', '#06b6d4']) # gray, green, red, cyan
    
    # Build RGBA array where the area outside the paddock (NaNs) is completely transparent (alpha=0)
    # This prevents the gray background from overlaying outside the paddock boundaries.
    tiers_gps_lookup = np.nan_to_num(tiers_gps, nan=0.0).astype(int)
    rgba_img = discrete_cmap(tiers_gps_lookup)
    rgba_img[np.isnan(tiers_gps), 3] = 0.0
    
    out_png_path = os.path.join(out_dir, "buffel_management_tiers.png")
    plt.imsave(out_png_path, rgba_img)
    print(f"[+] Saved colorized preview PNG to: {out_png_path}")
    
    # Calculate area statistics dynamically using pixel area in hectares
    transform = da_spike_rgb.rio.transform()
    pixel_width = abs(transform[0])
    pixel_height = abs(transform[4])
    pixel_area_ha = (pixel_width * pixel_height) / 10000.0
    
    area_t1 = np.sum(tiers == 1) * pixel_area_ha
    area_t2 = np.sum(tiers == 2) * pixel_area_ha
    area_t3 = np.sum(tiers == 3) * pixel_area_ha
    area_total = np.sum(valid_mask) * pixel_area_ha
    
    print("\n" + "=" * 40)
    print("         MANAGEMENT TIER AREAS")
    print("=" * 40)
    print(f"Total Analysis Area:  {area_total:.2f} ha")
    print(f"Tier 1 (Eradication): {area_t1:.2f} ha ({area_t1/area_total*100:.1f}%)")
    print(f"Tier 2 (Containment): {area_t2:.2f} ha ({area_t2/area_total*100:.1f}%)")
    print(f"Tier 3 (Refugia):     {area_t3:.2f} ha ({area_t3/area_total*100:.1f}%)")
    print("=" * 40)
    
    # Save stats to CSV
    stats_df = pd.DataFrame([{
        'Tier_1_Eradication_Ha': area_t1,
        'Tier_2_Containment_Ha': area_t2,
        'Tier_3_Refugia_Ha': area_t3,
        'Total_Area_Ha': area_total
    }])
    stats_csv = os.path.join(out_dir, "management_tier_statistics.csv")
    stats_df.to_csv(stats_csv, index=False)
    
    # 7. Generate Leaflet HTML Map with Tiers Layer Overlay
    bounds = da_tiers_gps.rio.bounds()
    leaflet_bounds = [[bounds[1], bounds[0]], [bounds[3], bounds[2]]]
    
    # Get centroid in EPSG:4326 for map center
    centroid = geom_multi.centroid
    center_lat = centroid.y
    center_lon = centroid.x
    
    # Add polygon outlines to Leaflet using GeoJSON directly
    geojson_data = gdf_wgs84.to_json()
    
    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <title>Keewong Property - 3-Tier Management Map</title>
    <meta charset="utf-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css" />
    <style>
        html, body {{
            height: 100%;
            margin: 0;
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
            background-color: #1a1a1a;
            color: #f3f3f3;
        }}
        #map {{
            width: 100%;
            height: 100%;
            position: absolute;
            top: 0;
            bottom: 0;
            z-index: 1;
        }}
        .info-panel {{
            padding: 12px;
            background: rgba(26, 26, 26, 0.95);
            box-shadow: 0 4px 15px rgba(0,0,0,0.5);
            border-radius: 8px;
            border: 1px solid #444;
            color: #f3f3f3;
            max-width: 300px;
        }}
        .info-panel h2 {{
            margin: 0 0 8px 0;
            font-size: 15px;
            color: #38bdf8;
            border-bottom: 1px solid #444;
            padding-bottom: 4px;
        }}
        .info-panel table {{
            width: 100%;
            font-size: 12px;
            border-collapse: collapse;
            margin-top: 6px;
        }}
        .info-panel td {{
            padding: 4px 0;
        }}
        .info-panel .val {{
            font-weight: bold;
            text-align: right;
        }}
        .legend {{
            padding: 8px;
            background: rgba(26, 26, 26, 0.95);
            border: 1px solid #444;
            border-radius: 8px;
            color: #f3f3f3;
            font-size: 11px;
            line-height: 16px;
        }}
        .legend i {{
            width: 14px;
            height: 14px;
            float: left;
            margin-right: 6px;
            border-radius: 2px;
            opacity: 0.85;
        }}
    </style>
</head>
<body>

<div id="map"></div>

<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
<script>
    var map = L.map('map').setView([{center_lat}, {center_lon}], 14);

    var satellite = L.tileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{{z}}/{{y}}/{{x}}', {{
        attribution: 'Tiles &copy; Esri'
    }}).addTo(map);

    var boundaryData = {geojson_data};
    
    // Draw paddock boundaries using standard Leaflet geoJSON
    L.geoJSON(boundaryData, {{
        style: function(feature) {{
            return {{
                color: '#ffffff',
                fillColor: '#000000',
                fillOpacity: 0.0,
                weight: 2,
                dashArray: '4, 4'
            }};
        }}
    }}).addTo(map);

    // Add Tier Layer Overlay
    var tierOverlay = L.imageOverlay('buffel_management_tiers.png', {leaflet_bounds}, {{
        opacity: 0.70,
        interactive: false
    }}).addTo(map);

    // Add legend
    var legend = L.control({{ position: 'bottomright' }});
    legend.onAdd = function (map) {{
        var div = L.DomUtil.create('div', 'legend');
        div.innerHTML = '<h3 style="margin:0 0 5px 0; font-size:12px; color:#38bdf8;">Weed Management Tiers</h3>' +
            '<i style="background: #ef4444"></i> <strong>Tier 2: Containment</strong> (Heavy Infestations)<br>' +
            '<i style="background: #22c55e"></i> <strong>Tier 1: Eradication</strong> (Light Infestations)<br>' +
            '<i style="background: #06b6d4"></i> <strong>Tier 3: Refugia</strong> (Stable Conservation Trees)<br>' +
            '<i style="background: #4b5563"></i> <strong>Background</strong> (Bare Soil / Native Grasses)<br>';
        return div;
    }};
    legend.addTo(map);

    // Add info control with exact hectares
    var info = L.control({{ position: 'bottomleft' }});
    info.onAdd = function(map) {{
        var div = L.DomUtil.create('div', 'info-panel');
        div.innerHTML = '<h2>3-Tier Weed Management Overlay</h2>' +
            '<table>' +
            '<tr><td>Total Buffered Area:</td><td class="val">{area_total:.1f} ha</td></tr>' +
            '<tr style="color: #22c55e"><td>Tier 1 (Eradication):</td><td class="val">{area_t1:.1f} ha ({area_t1/area_total*100:.1f}%)</td></tr>' +
            '<tr style="color: #ef4444"><td>Tier 2 (Containment):</td><td class="val">{area_t2:.1f} ha ({area_t2/area_total*100:.1f}%)</td></tr>' +
            '<tr style="color: #06b6d4"><td>Tier 3 (Refugia):</td><td class="val">{area_t3:.1f} ha ({area_t3/area_total*100:.1f}%)</td></tr>' +
            '</table>' +
            '<p style="font-size:10px; margin-top:8px; color:#888; border-top:1px solid #444; padding-top:4px;">' +
            'This layer is dynamically generated using MSAVI standard deviation (annual volatility) and targeted rain-event response signatures.' +
            '</p>';
        return div;
    }};
    info.addTo(map);

</script>
</body>
</html>
"""
    
    html_path = os.path.join(out_dir, "keewong_management_tiers_map.html")
    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    print(f"[+] Saved management tiers Leaflet HTML map to: {html_path}")
    
    print("\n" + "=" * 60)
    print("[+] 3-TIER WEED MANAGEMENT LAYER GENERATION COMPLETE!")
    print("=" * 60)

if __name__ == "__main__":
    main()
