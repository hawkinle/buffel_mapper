import os
# Set GDAL network timeouts to prevent hanging on remote HTTP downloads
os.environ["GDAL_HTTP_TIMEOUT"] = "15"
os.environ["GDAL_HTTP_CONNECTTIMEOUT"] = "5"
os.environ["GDAL_HTTP_MAX_RETRY"] = "3"

import requests
# Patch requests to have a default timeout of 15 seconds to prevent hanging on STAC API / signing calls
original_request = requests.Session.request
def timeout_request(self, method, url, **kwargs):
    if 'timeout' not in kwargs or kwargs['timeout'] is None:
        kwargs['timeout'] = 15
    return original_request(self, method, url, **kwargs)
requests.Session.request = timeout_request

import argparse
import datetime
import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd
import shapely.geometry
from shapely.geometry import box
import rioxarray
import pystac_client
import planetary_computer
import rasterio
import xml.etree.ElementTree as ET
import xarray as xr
import pandas as pd

def load_boundary_geometry(file_path):
    import geopandas as gpd
    import shapely.geometry
    
    ext = os.path.splitext(file_path)[1].lower()
    
    if ext == ".kml":
        print(f"[*] Parsing KML boundary from: {file_path}")
        tree = ET.parse(file_path)
        root = tree.getroot()
        namespaces = {
            'kml': 'http://www.opengis.net/kml/2.2',
            'gx': 'http://www.google.com/kml/ext/2.2'
        }
        coordinates_elements = root.findall('.//kml:coordinates', namespaces)
        if not coordinates_elements:
            coordinates_elements = root.findall('.//coordinates')
        if not coordinates_elements:
            raise ValueError("Could not find coordinates in KML.")
        coord_text = coordinates_elements[0].text.strip()
        coords = []
        for pt in coord_text.split():
            parts = pt.split(',')
            if len(parts) >= 2:
                lon, lat = float(parts[0]), float(parts[1])
                coords.append((lon, lat))
        return shapely.geometry.Polygon(coords)
        
    elif ext in [".shp", ".geojson", ".gpkg"]:
        print(f"[*] Reading vector boundary shapefile from: {file_path}")
        gdf = gpd.read_file(file_path)
        if gdf.empty:
            raise ValueError("Vector file is empty.")
        gdf_wgs84 = gdf.to_crs("EPSG:4326")
        geom = gdf_wgs84.geometry.unary_union
        return geom
    else:
        raise ValueError(f"Unsupported boundary file format: {ext}")

def get_leaflet_coordinates(polygon):
    import shapely.geometry
    if isinstance(polygon, shapely.geometry.MultiPolygon):
        coords = []
        for poly in polygon.geoms:
            poly_coords = []
            x_coords, y_coords = poly.exterior.coords.xy
            for x, y in zip(x_coords, y_coords):
                poly_coords.append([y, x])
            coords.append(poly_coords)
        return coords
    else:
        x_coords, y_coords = polygon.exterior.coords.xy
        coords = []
        for x, y in zip(x_coords, y_coords):
            coords.append([y, x]) # Leaflet expects [latitude, longitude]
        return coords

def get_clip_geometry(polygon, raster_crs, buffer_pixels=0):
    import geopandas as gpd
    gdf_poly = gpd.GeoDataFrame(geometry=[polygon], crs="EPSG:4326")
    gdf_poly_proj = gdf_poly.to_crs(raster_crs)
    
    if buffer_pixels > 0:
        buffer_meters = buffer_pixels * 10.0
        # Buffer in projected coordinates (meters)
        buffered_geom = gdf_poly_proj.geometry.buffer(buffer_meters).iloc[0]
        return buffered_geom, raster_crs
    else:
        return polygon, "EPSG:4326"

def parse_rain_file(rain_path, start_year=2016, warm_season_only=True, cluster_days=14):
    print(f"[*] Parsing rain dates from: {rain_path}")
    raw_events = []
    with open(rain_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or 'Dates with' in line:
                continue
            parts = line.split(' - ')
            if len(parts) == 2:
                try:
                    d = datetime.datetime.strptime(parts[0], '%d/%m/%Y')
                    if d.year >= start_year:
                        if warm_season_only and d.month not in [10, 11, 12, 1, 2, 3, 4]:
                            continue
                        amount_str = parts[1].replace(' mm', '').strip()
                        amount = float(amount_str)
                        raw_events.append((d.date(), amount))
                except ValueError:
                    continue
    raw_events = sorted(raw_events, key=lambda x: x[0])
    season_label = "warm-season (Oct-Apr)" if warm_season_only else "all-season"
    print(f"[+] Found {len(raw_events)} raw {season_label} rain records of >= 15mm since {start_year}.")

    if not raw_events or cluster_days <= 0:
        return [(d, d, amt) for d, amt in raw_events]

    # Group consecutive rain dates within `cluster_days` of cluster start into a single pulse
    clustered_events = []
    cluster_start, cluster_end, cluster_sum = raw_events[0][0], raw_events[0][0], raw_events[0][1]
    for d, amt in raw_events[1:]:
        if (d - cluster_start).days <= cluster_days:
            cluster_end = d
            cluster_sum += amt
        else:
            clustered_events.append((cluster_start, cluster_end, round(cluster_sum, 1)))
            cluster_start, cluster_end, cluster_sum = d, d, amt
    clustered_events.append((cluster_start, cluster_end, round(cluster_sum, 1)))
    print(f"[+] Grouped into {len(clustered_events)} distinct rain-pulse clusters ({cluster_days}-day window).")
    return clustered_events

def dilate_mask_2d(mask_2d, radius=2):
    """Dilate a 2D boolean mask by `radius` pixels using pure NumPy to buffer cloud/shadow edges."""
    if radius <= 0:
        return mask_2d.astype(bool)
    padded = np.pad(mask_2d.astype(bool), radius, mode='edge')
    h, w = mask_2d.shape
    out = np.zeros((h, w), dtype=bool)
    for dy in range(-radius, radius + 1):
        for dx in range(-radius, radius + 1):
            out |= padded[radius + dy : radius + dy + h, radius + dx : radius + dx + w]
    return out

def save_raw_float_tiff(data_2d, out_path, crs, transform):
    """Save a 2D float array as a lossless single-band float32 GeoTIFF."""
    meta = {
        'driver': 'GTiff',
        'dtype': 'float32',
        'nodata': np.nan,
        'width': data_2d.shape[1],
        'height': data_2d.shape[0],
        'count': 1,
        'crs': crs,
        'transform': transform
    }
    with rasterio.open(out_path, 'w', **meta) as dst:
        dst.write(data_2d.astype(np.float32), 1)

def save_transparent_png(data_2d, vmin, vmax, cmap_name, out_path):
    """Save a 2D float array as a colorized RGBA PNG where NaN pixels are completely transparent."""
    cmap = plt.get_cmap(cmap_name)
    valid = ~np.isnan(data_2d)
    norm = np.zeros_like(data_2d, dtype=np.float32)
    norm[valid] = np.clip((data_2d[valid] - vmin) / (vmax - vmin), 0.0, 1.0)
    rgba = cmap(norm)
    rgba[~valid, 3] = 0.0
    plt.imsave(out_path, rgba)

def calculate_index(red, nir, mask, index_name="ndvi"):
    r = red.where(mask)
    n = nir.where(mask)
    if index_name.lower() == "msavi":
        # MSAVI formula: (2 * NIR + 1 - sqrt((2 * NIR + 1)^2 - 8 * (NIR - Red))) / 2
        sqrt_term = np.clip((2 * n + 1)**2 - 8 * (n - r), 0, None)
        return (2 * n + 1 - np.sqrt(sqrt_term)) / 2
    else:
        # Default to NDVI
        return (n - r) / (n + r)

def run_targeted_analysis(kml_path, rain_path, out_dir, threshold=0.15, max_cloud=10.0, index_name="msavi", buffer_pixels=0, warm_season_only=True):
    os.makedirs(out_dir, exist_ok=True)
    polygon = load_boundary_geometry(kml_path)
    bbox_coords = list(polygon.bounds)
    
    if buffer_pixels > 0:
        # Buffer the search bounding box by ~30m per pixel in degrees to ensure all scenes are queried
        degree_buffer = buffer_pixels * 0.0003
        bbox_coords = [
            bbox_coords[0] - degree_buffer,
            bbox_coords[1] - degree_buffer,
            bbox_coords[2] + degree_buffer,
            bbox_coords[3] + degree_buffer
        ]
        
    rain_events = parse_rain_file(rain_path, warm_season_only=warm_season_only)
    
    # 1. Connect to Planetary Computer
    print("[*] Connecting to Planetary Computer STAC API...")
    catalog = pystac_client.Client.open(
        "https://planetarycomputer.microsoft.com/api/stac/v1",
        modifier=planetary_computer.sign_inplace
    )
    
    # 2. Search all scenes for the decade to index them in memory (optimizes query time)
    start_date = "2016-01-01"
    end_date = datetime.date.today().strftime("%Y-%m-%d")
    print(f"[*] Querying full Sentinel-2 catalog for the property ({start_date} to {end_date})...")
    
    search = catalog.search(
        collections=["sentinel-2-l2a"],
        bbox=bbox_coords,
        datetime=f"{start_date}/{end_date}",
        query={"eo:cloud_cover": {"lt": max_cloud}}
    )
    
    items = search.item_collection()
    if len(items) == 0:
        print("[!] No scenes found.")
        return
        
    print(f"[+] Loaded {len(items)} scenes into memory. Indexing by date...")
    
    # Index scenes by date
    indexed_scenes = []
    for item in items:
        dt = datetime.datetime.strptime(item.properties["datetime"][:10], "%Y-%m-%d").date()
        indexed_scenes.append((dt, item))
    # Sort chronologically
    indexed_scenes = sorted(indexed_scenes, key=lambda x: x[0])
    
    geom = box(*bbox_coords)
    gdf = gpd.GeoDataFrame(geometry=[geom], crs="EPSG:4326")
    
    diff_list = []
    pre_list = []
    event_details = []
    
    crs_ref = None
    transform_ref = None
    leaflet_bounds = None
    ref_template = None
    
    # 3. Analyze each rain event cluster
    print("\n[*] Running targeted pre/post rain analysis for each event pulse...")
    for idx, (rain_start, rain_end, amount) in enumerate(rain_events):
        pulse_str = f"{rain_start.strftime('%Y-%m-%d')}" if rain_start == rain_end else f"{rain_start.strftime('%Y-%m-%d')} to {rain_end.strftime('%Y-%m-%d')}"
        print(f"    [{idx+1}/{len(rain_events)}] Pulse: {pulse_str} ({amount} mm)...")
        
        # Pre-rain window: 1 to 12 days before the start of the rain pulse
        pre_min = rain_start - datetime.timedelta(days=12)
        pre_max = rain_start - datetime.timedelta(days=1)
        
        # Post-rain window: 12 to 28 days after the rain pulse start (captures C4 green-up peak)
        post_min = rain_start + datetime.timedelta(days=12)
        post_max = rain_end + datetime.timedelta(days=28)
        
        # Find best pre-rain scene
        pre_scenes = [s for s in indexed_scenes if pre_min <= s[0] <= pre_max]
        if not pre_scenes:
            # Try slightly wider window 1-15 days
            pre_min_wide = rain_start - datetime.timedelta(days=15)
            pre_scenes = [s for s in indexed_scenes if pre_min_wide <= s[0] <= pre_max]
            
        # Find best post-rain scene
        post_scenes = [s for s in indexed_scenes if post_min <= s[0] <= post_max]
        if not post_scenes:
            # Try slightly wider window 10-35 days
            post_max_wide = rain_end + datetime.timedelta(days=35)
            post_scenes = [s for s in indexed_scenes if post_min <= s[0] <= post_max_wide]
            
        if not pre_scenes or not post_scenes:
            print("      [-] Skipped: No matching pre-rain or post-rain scenes found in STAC cache.")
            continue
            
        # Select least cloudy scenes in the windows
        pre_item = min(pre_scenes, key=lambda x: x[1].properties["eo:cloud_cover"])[1]
        post_item = min(post_scenes, key=lambda x: x[1].properties["eo:cloud_cover"])[1]
        
        print(f"      [+] Selected Pre-Rain: {pre_item.properties['datetime'][:10]} (Cloud: {pre_item.properties['eo:cloud_cover']:.1f}%)")
        print(f"      [+] Selected Post-Rain: {post_item.properties['datetime'][:10]} (Cloud: {post_item.properties['eo:cloud_cover']:.1f}%)")
        
        try:
            # Calculate Pre-Rain Index
            pre_red_url = pre_item.assets["B04"].href
            pre_nir_url = pre_item.assets["B08"].href
            pre_scl_url = pre_item.assets["SCL"].href
            
            da_red_pre = rioxarray.open_rasterio(pre_red_url)
            gdf_proj = gdf.to_crs(da_red_pre.rio.crs)
            proj_bbox = gdf_proj.total_bounds
            
            clip_geom, clip_crs = get_clip_geometry(polygon, da_red_pre.rio.crs, buffer_pixels)
            red_pre_crop = da_red_pre.rio.clip_box(*proj_bbox).rio.clip([clip_geom], crs=clip_crs, all_touched=True)
            red_pre = red_pre_crop.where(red_pre_crop != da_red_pre.rio.nodata).astype(np.float32) / 10000.0
            
            da_nir_pre = rioxarray.open_rasterio(pre_nir_url)
            nir_pre_crop = da_nir_pre.rio.clip_box(*proj_bbox).rio.clip([clip_geom], crs=clip_crs, all_touched=True)
            nir_pre = nir_pre_crop.where(nir_pre_crop != da_nir_pre.rio.nodata).astype(np.float32) / 10000.0
            
            da_scl_pre = rioxarray.open_rasterio(pre_scl_url)
            scl_pre_crop = da_scl_pre.rio.clip_box(*proj_bbox).rio.clip([clip_geom], crs=clip_crs, all_touched=True)
            scl_pre_res = scl_pre_crop.rio.reproject_match(red_pre, resampling=rasterio.enums.Resampling.nearest)
            scl_pre_val = scl_pre_res.squeeze()
            
            invalid_scl = [0, 1, 3, 8, 9, 10, 11]
            cloud_pre = dilate_mask_2d(np.isin(scl_pre_val.values, invalid_scl), radius=2)
            valid_mask_pre = ~cloud_pre
            
            index_pre = calculate_index(red_pre, nir_pre, valid_mask_pre, index_name)
            index_pre = index_pre.squeeze()
            
            # Calculate Post-Rain Index
            post_red_url = post_item.assets["B04"].href
            post_nir_url = post_item.assets["B08"].href
            post_scl_url = post_item.assets["SCL"].href
            
            da_red_post = rioxarray.open_rasterio(post_red_url)
            clip_geom_post, clip_crs_post = get_clip_geometry(polygon, da_red_post.rio.crs, buffer_pixels)
            red_post_crop = da_red_post.rio.clip_box(*proj_bbox).rio.clip([clip_geom_post], crs=clip_crs_post, all_touched=True)
            red_post = red_post_crop.where(red_post_crop != da_red_post.rio.nodata).astype(np.float32) / 10000.0
            
            da_nir_post = rioxarray.open_rasterio(post_nir_url)
            nir_post_crop = da_nir_post.rio.clip_box(*proj_bbox).rio.clip([clip_geom_post], crs=clip_crs_post, all_touched=True)
            nir_post = nir_post_crop.where(nir_post_crop != da_nir_post.rio.nodata).astype(np.float32) / 10000.0
            
            da_scl_post = rioxarray.open_rasterio(post_scl_url)
            scl_post_crop = da_scl_post.rio.clip_box(*proj_bbox).rio.clip([clip_geom_post], crs=clip_crs_post, all_touched=True)
            scl_post_res = scl_post_crop.rio.reproject_match(red_post, resampling=rasterio.enums.Resampling.nearest)
            scl_post_val = scl_post_res.squeeze()
            
            cloud_post = dilate_mask_2d(np.isin(scl_post_val.values, invalid_scl), radius=2)
            valid_mask_post = ~cloud_post
            
            index_post = calculate_index(red_post, nir_post, valid_mask_post, index_name)
            index_post = index_post.squeeze()
            
            # Calculate growth difference
            index_diff = index_post - index_pre
            
            # Align to reference template if needed
            if ref_template is None:
                ref_template = index_diff
                crs_ref = index_diff.rio.crs
                transform_ref = index_diff.rio.transform()
            else:
                index_diff = index_diff.rio.reproject_match(ref_template)
                index_pre = index_pre.rio.reproject_match(ref_template)
                
            diff_list.append(index_diff)
            pre_list.append(index_pre)
                
            event_details.append({
                'Rain_Date': rain_start,
                'Rain_End_Date': rain_end,
                'Pre_Date': pre_item.properties['datetime'][:10],
                'Post_Date': post_item.properties['datetime'][:10],
                'Rain_Amt_mm': amount
            })
            print("      [+] Successfully calculated growth difference.")
            
        except Exception as e:
            print(f"      [!] Failed to process scene files: {e}")
            continue
            
    if len(diff_list) == 0:
        print("[!] No rain events successfully analyzed.")
        return
        
    print(f"\n[+] Successfully analyzed {len(diff_list)} rain pulse events over the decade.")
    
    # 4. Stack and compute statistics
    print("[*] Stacking events and calculating targeted metrics...")
    diff_stack = xr.concat(diff_list, dim='event')
    pre_stack = xr.concat(pre_list, dim='event')
    
    diff_vals = diff_stack.values  # (n_events, h, w)
    pre_vals = pre_stack.values    # (n_events, h, w)
    
    # True footprint mask inside the buffered boundary
    footprint_mask = np.any(~np.isnan(diff_vals), axis=0)
    
    # Metric 1: 90th-Percentile Growth Spike (robust against single-scene cloud halo artifacts)
    with np.errstate(all='ignore'):
        max_spike_2d = np.nanpercentile(diff_vals, 90, axis=0).astype(np.float32)
    max_spike_2d[~footprint_mask] = np.nan
    
    # Metric 2: True Growth Response Frequency
    # Count only events where the pixel was observed AND not already saturated green prior to rain (pre_vals < 0.25)
    valid_dry_events = (~np.isnan(diff_vals)) & (pre_vals < 0.25)
    spiked_events = valid_dry_events & (diff_vals >= threshold)
    
    valid_event_count = np.sum(valid_dry_events, axis=0)
    spiked_count = np.sum(spiked_events, axis=0)
    
    with np.errstate(divide='ignore', invalid='ignore'):
        growth_frequency_2d = np.where(
            valid_event_count > 0,
            spiked_count / valid_event_count,
            0.0
        ).astype(np.float32)
    growth_frequency_2d[~footprint_mask] = np.nan
    
    # 5. Save lossless single-band float32 GeoTIFFs
    max_tif_path = os.path.join(out_dir, "buffel_max_growth_spike.tif")
    freq_tif_path = os.path.join(out_dir, "buffel_growth_frequency.tif")
    
    save_raw_float_tiff(max_spike_2d, max_tif_path, crs_ref, transform_ref)
    save_raw_float_tiff(growth_frequency_2d, freq_tif_path, crs_ref, transform_ref)
    
    # Reproject to EPSG:4326 for accurate transparent web map PNG overlays
    da_spike_out = xr.DataArray(max_spike_2d, coords={'y': ref_template.y, 'x': ref_template.x}, dims=['y', 'x']).rio.write_crs(crs_ref)
    da_freq_out = xr.DataArray(growth_frequency_2d, coords={'y': ref_template.y, 'x': ref_template.x}, dims=['y', 'x']).rio.write_crs(crs_ref)
    
    da_spike_gps = da_spike_out.rio.reproject("EPSG:4326")
    da_freq_gps = da_freq_out.rio.reproject("EPSG:4326")
    bounds = da_spike_gps.rio.bounds()
    leaflet_bounds = [[bounds[1], bounds[0]], [bounds[3], bounds[2]]]
    
    max_png_path = os.path.join(out_dir, "buffel_max_growth_spike.png")
    freq_png_path = os.path.join(out_dir, "buffel_growth_frequency.png")
    spike_vmax = 0.20 if index_name.lower() == "msavi" else 0.35
    freq_vmax = 0.25 if index_name.lower() == "msavi" else 0.50
    save_transparent_png(da_spike_gps.values.squeeze(), 0.0, spike_vmax, "YlGn", max_png_path)
    save_transparent_png(da_freq_gps.values.squeeze(), 0.0, freq_vmax, "plasma", freq_png_path)
    
    # Save report plot
    report_path = os.path.join(out_dir, "targeted_analysis_report.png")
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    
    im0 = axes[0].imshow(max_spike_2d, cmap="YlGn", vmin=0.0, vmax=spike_vmax)
    axes[0].set_title(f"90th-Percentile Growth Spike ({index_name.upper()} Diff)\n(Dark green = Vigorous warm-season rain response)")
    axes[0].axis("off")
    fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)
    
    im1 = axes[1].imshow(growth_frequency_2d, cmap="plasma", vmin=0.0, vmax=freq_vmax)
    axes[1].set_title(f"Growth Response Frequency ({index_name.upper()} Diff >= {threshold})\n(Yellow/bright = Consistently spikes after warm rain)")
    axes[1].axis("off")
    fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)
    
    plt.tight_layout()
    plt.savefig(report_path, dpi=150)
    plt.close()
    
    # 6. Save event details to CSV
    df_events = pd.DataFrame(event_details)
    csv_path = os.path.join(out_dir, "analyzed_rain_events.csv")
    df_events.to_csv(csv_path, index=False)
    print(f"[+] Saved analyzed event log to: {csv_path}")
    
    # 7. Generate Leaflet HTML Map
    leaflet_coords = get_leaflet_coordinates(polygon)
    
    # Flatten coordinates to compute center
    flat_lats = []
    flat_lons = []
    if isinstance(leaflet_coords[0][0], list):
        for poly_coords in leaflet_coords:
            for c in poly_coords:
                flat_lats.append(c[0])
                flat_lons.append(c[1])
    else:
        for c in leaflet_coords:
            flat_lats.append(c[0])
            flat_lons.append(c[1])
            
    center_lat = sum(flat_lats) / len(flat_lats)
    center_lon = sum(flat_lons) / len(flat_lons)
    
    # Compute property area stats
    # 1 pixel = 0.01 ha
    valid_mask = ~np.isnan(max_spike_2d)
    total_area_ha = np.sum(valid_mask) * 0.01
    
    # High-confidence core Buffel grass (calibrated for 90th-percentile MSAVI vs NDVI)
    spike_core_thresh = 0.105 if index_name.lower() == "msavi" else 0.20
    freq_core_thresh = 0.05 if index_name.lower() == "msavi" else 0.20
    high_confidence_mask = (max_spike_2d >= spike_core_thresh) & (growth_frequency_2d >= freq_core_thresh) & valid_mask
    high_conf_area_ha = np.sum(high_confidence_mask) * 0.01
    infest_pct = (high_conf_area_ha / total_area_ha) * 100 if total_area_ha > 0 else 0
    
    print(f"[+] High-confidence Buffel infestation area: {high_conf_area_ha:.2f} ha ({infest_pct:.1f}%)")
    
    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <title>Keewong Property - Targeted Rain-Event Weed Map</title>
    <meta charset="utf-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css" />
    <style>
        html, body {{
            height: 100%;
            margin: 0;
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
            background-color: #1e1e24;
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
            padding: 15px;
            background: rgba(30, 30, 36, 0.95);
            box-shadow: 0 4px 15px rgba(0,0,0,0.5);
            border-radius: 8px;
            border: 1px solid #444;
            color: #f3f3f3;
            max-width: 300px;
        }}
        .info-panel h2 {{
            margin: 0 0 10px 0;
            font-size: 15px;
            color: #4ade80;
            border-bottom: 1px solid #444;
            padding-bottom: 5px;
        }}
        .stat-row {{
            display: flex;
            justify-content: space-between;
            margin: 6px 0;
            font-size: 13px;
        }}
        .stat-val {{
            font-weight: bold;
            color: #e2e8f0;
        }}
        .stat-highlight {{
            color: #f87171;
            font-weight: bold;
        }}
        .legend {{
            padding: 10px;
            background: rgba(30, 30, 36, 0.95);
            border: 1px solid #444;
            border-radius: 8px;
            line-height: 18px;
            color: #f3f3f3;
            font-size: 12px;
            max-width: 200px;
        }}
        .legend i {{
            width: 18px;
            height: 18px;
            float: left;
            margin-right: 8px;
            opacity: 0.85;
            border-radius: 2px;
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

    var boundaryCoords = {leaflet_coords};
    var boundary = L.polygon(boundaryCoords, {{
        color: '#38bdf8',
        fillColor: 'transparent',
        weight: 3,
        opacity: 0.9
    }}).addTo(map);

    // Overlay 1: Maximum Growth Spike (YlGn colormap)
    var maxSpikeOverlay = L.imageOverlay('buffel_max_growth_spike.png', {leaflet_bounds}, {{
        opacity: 0.65,
        interactive: false
    }}).addTo(map);

    // Overlay 2: Growth Response Frequency (Plasma colormap)
    var freqOverlay = L.imageOverlay('buffel_growth_frequency.png', {leaflet_bounds}, {{
        opacity: 0.65,
        interactive: false
    }});

    var basemaps = {{
        "Satellite Imagery": satellite
    }};
    
    var overlayMaps = {{
        "Property Boundary": boundary,
        "Max Growth Spike (Greenness)": maxSpikeOverlay,
        "Growth Response Frequency (Plasma)": freqOverlay
    }};
    
    L.control.layers(basemaps, overlayMaps, {{ collapsed: false, position: 'topright' }}).addTo(map);

    // Info panel
    var info = L.control({{ position: 'bottomleft' }});
    info.onAdd = function (map) {{
        var div = L.DomUtil.create('div', 'info-panel');
        div.innerHTML = '<h2>Rain-Event Targeted Analysis</h2>' +
            '<div class="stat-row"><span>Total Property:</span><span class="stat-val">{total_area_ha:.1f} ha</span></div>' +
            '<div class="stat-row"><span>Analyzed Rain Events:</span><span class="stat-val">{len(diff_list)} events</span></div>' +
            '<div class="stat-row" style="border-top:1px solid #444; margin-top:8px; padding-top:8px;">' +
            '<span>High-Confidence Buffel:</span><span class="stat-highlight">{high_conf_area_ha:.1f} ha</span></div>' +
            '<div class="stat-row"><span>Infestation Rate:</span><span class="stat-highlight">{infest_pct:.1f}%</span></div>' +
            '<div class="stat-row" style="font-size: 11px; color:#aaa; margin-top:10px; border-top:1px solid #444; padding-top:5px;">' +
            '<span>Method:</span><span>{index_name.upper()} Spike Pre- vs Post-Rain</span></div>';
        return div;
    }};
    info.addTo(map);

    // Legend
    var legend = L.control({{ position: 'bottomright' }});
    legend.onAdd = function (map) {{
        var div = L.DomUtil.create('div', 'legend');
        div.innerHTML = '<h3 style="margin:0 0 5px 0; font-size:12px; color:#4ade80;">Max Growth Spike</h3>' +
            '<i style="background: #00441b"></i> Strong Growth Response (High Prob)<br>' +
            '<i style="background: #a1d99b"></i> Moderate Response<br>' +
            '<i style="background: #f7fcf5"></i> No Growth (Soil / Stable Tree)<br>';
        return div;
    }};
    legend.addTo(map);

    map.fitBounds(boundary.getBounds());
</script>

</body>
</html>
"""
    
    html_path = os.path.join(out_dir, "keewong_targeted_weed_map.html")
    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    print(f"[+] Saved targeted Leaflet HTML map: {html_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Targeted Rain-Event Weed Analyzer")
    parser.add_argument("--kml", type=str, required=True, help="Path to KML or Shapefile boundary")
    parser.add_argument("--rain-file", type=str, required=True, help="Path to BoM rain text file")
    parser.add_argument("--out-dir", type=str, default="output_targeted", help="Output directory")
    parser.add_argument("--threshold", type=float, default=0.15, help="Index change threshold")
    parser.add_argument("--max-cloud", type=float, default=10.0, help="Max cloud cover %")
    parser.add_argument("--index", type=str, default="msavi", choices=["ndvi", "msavi"], help="Vegetation index to use (ndvi or msavi)")
    parser.add_argument("--buffer-pixels", type=int, default=0, help="Buffer boundaries by N pixels (1 pixel = 10m)")
    parser.add_argument("--include-winter", action="store_true", help="Include cool-season winter rain events (May-Sep) in addition to warm season (Oct-Apr)")
    
    args = parser.parse_args()
    run_targeted_analysis(
        args.kml,
        args.rain_file,
        args.out_dir,
        args.threshold,
        args.max_cloud,
        args.index,
        args.buffer_pixels,
        warm_season_only=not args.include_winter
    )
