#!/usr/bin/env python3
"""
Unified 3-Tier Weed Management Spatial Pipeline
Author: Antigravity Spatial Remote Sensing Pipeline
"""

import os
import argparse
import time
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
import rasterio
from rasterio.enums import Resampling
import rioxarray
import pystac_client
from pyproj import Transformer

def parse_args():
    parser = argparse.ArgumentParser(description="Unified 3-Tier Weed Management Spatial Pipeline")
    parser.add_argument("-s", "--shapefile", required=True, help="Path to boundary shapefile (.shp)")
    parser.add_argument("-k", "--max-spike", required=True, help="Path to targeted max growth spike raster (.tif)")
    parser.add_argument("-f", "--growth-freq", required=True, help="Path to growth frequency raster (.tif)")
    parser.add_argument("-o", "--output-dir", required=True, help="Directory to save final classified layers")
    parser.add_argument("--buffer-pixels", type=int, default=2, help="Buffer pixels around paddock boundaries (default 2)")
    parser.add_argument("--start-year", type=int, default=2023, help="Start year for dry-season baseline (default 2023)")
    parser.add_argument("--end-year", type=int, default=2025, help="End year for dry-season baseline (default 2025)")
    parser.add_argument("--cmap-name", default='viridis', help="Colormap name if input rasters are RGB colorized (default 'viridis')")
    parser.add_argument("--cmap-vmin", type=float, default=-0.1, help="Vmin if input rasters are RGB colorized (default -0.1)")
    parser.add_argument("--cmap-vmax", type=float, default=0.4, help="Vmax if input rasters are RGB colorized (default 0.4)")
    return parser.parse_args()

def read_raster_raw(path, cmap_name='viridis', vmin=-0.1, vmax=0.4):
    """
    Reads a raster. If it has 3+ bands, decodes it from RGB back to raw float values using a colormap lookup.
    If it is a single band, reads it directly as raw float values.
    """
    with rasterio.open(path) as src:
        bands = src.count
        data = src.read()
        transform = src.transform
        crs = src.crs
        
    if bands >= 3:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        
        print(f"    [*] Decoding RGB colorized raster back to float using colormap '{cmap_name}'...")
        cmap = plt.get_cmap(cmap_name)
        lut = np.array([cmap(i)[:3] for i in range(256)])
        
        h, w = data.shape[1], data.shape[2]
        flat_rgb = data[:3, :, :].reshape(3, -1).T
        flat_rgb = np.clip(flat_rgb, 0.0, 1.0)
        
        # Vectorized closest color matching
        dists = np.sum((flat_rgb[:, None, :] - lut[None, :, :])**2, axis=2)
        lut_indices = np.argmin(dists, axis=1)
        normalized = lut_indices.reshape(h, w).astype(float) / 255.0
        return vmin + normalized * (vmax - vmin), transform, crs
    else:
        print("    [*] Reading single-band raw float raster directly...")
        return data[0].astype(float), transform, crs

def generate_leaflet_map(tiers_wgs84, gdf_wgs84, output_html, bounds):
    """
    Generates a standalone Leaflet interactive map with paddock boundaries and colorized transparent tiers layer.
    """
    import json
    
    # Define color scheme (1: Green, 2: Red, 3: Cyan)
    colormap_js = """
    function getColor(val) {
        if (val === 1) return '#2ca02c'; // Tier 1 - Eradication (Green)
        if (val === 2) return '#d62728'; // Tier 2 - Containment (Red)
        if (val === 3) return '#00ffff'; // Tier 3 - Refugia (Cyan)
        return 'transparent';
    }
    """
    
    # Calculate bounding box zoom center
    center_lat = (bounds[1] + bounds[3]) / 2.0
    center_lon = (bounds[0] + bounds[2]) / 2.0
    
    # Extract geometries as GeoJSON
    geojson_str = gdf_wgs84.to_json()
    
    # Convert classified tiers to a tiny GeoJSON grid or render as image overlay
    # For high performance and standalone HTML, we will serialize the raster cells to a Leaflet layer
    # or generate a direct Leaflet map structure.
    
    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <title>3-Tier Weed Management Map</title>
    <meta charset="utf-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css" />
    <script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
    <style>
        html, body, #map {{
            width: 100%;
            height: 100%;
            margin: 0;
            padding: 0;
            background: #111;
        }}
        .legend {{
            padding: 10px;
            font: 14px/16px Arial, Helvetica, sans-serif;
            background: rgba(30, 30, 30, 0.9);
            color: #fff;
            box-shadow: 0 0 15px rgba(0,0,0,0.5);
            border-radius: 5px;
            line-height: 20px;
            border: 1px solid #444;
        }}
        .legend i {{
            width: 18px;
            height: 18px;
            float: left;
            margin-right: 8px;
            opacity: 0.8;
            border-radius: 3px;
        }}
    </style>
</head>
<body>
    <div id="map"></div>
    <script>
        var map = L.map('map').setView([{center_lat}, {center_lon}], 13);
        
        // Dark theme basemap
        L.tileLayer('https://{{s}}.basemaps.cartocdn.com/dark_all/{{z}}/{{x}}/{{y}}{{r}}.png', {{
            attribution: '&copy; OpenStreetMap &copy; CARTO',
            subdomains: 'abcd',
            maxZoom: 20
        }}).addTo(map);
        
        // Add paddock boundaries
        var paddockGeoJSON = {geojson_str};
        L.geoJSON(paddockGeoJSON, {{
            style: {{
                color: "#ff7800",
                weight: 2,
                opacity: 0.8,
                fillOpacity: 0.05
            }}
        }}).addTo(map);
        
        // Add classified layers as transparent image overlay
        // Image overlay is generated directly by the python script as a transparent WGS84 PNG
        var imageUrl = 'buffel_management_tiers.png';
        var imageBounds = [[{bounds[1]}, {bounds[0]}], [{bounds[3]}, {bounds[2]}]];
        L.imageOverlay(imageUrl, imageBounds, {{
            opacity: 0.85,
            interactive: false
        }}).addTo(map);
        
        // Legend
        var legend = L.control({{position: 'bottomright'}});
        legend.onAdd = function (map) {{
            var div = L.DomUtil.create('div', 'legend'),
                grades = [1, 2, 3],
                labels = [
                    "Tier 1 - Eradication (Light/Outlier)", 
                    "Tier 2 - Containment (Core Heavy)", 
                    "Tier 3 - Native Refugia (Evergreen Tree)"
                ],
                colors = ["#2ca02c", "#d62728", "#00ffff"];

            div.innerHTML = '<strong>3-Tier Management Layer</strong><br><br>';
            for (var i = 0; i < grades.length; i++) {{
                div.innerHTML +=
                    '<i style="background:' + colors[i] + '"></i> ' + labels[i] + '<br>';
            }}
            return div;
        }};
        legend.addTo(map);
    </script>
</body>
</html>
"""
    with open(output_html, 'w') as f:
        f.write(html_content)

def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    
    print("=============================================================")
    print("[*] UNIFIED 3-TIER WEED MANAGEMENT PIPELINE STARTING")
    print("=============================================================")
    
    # 1. Load Paddock Boundaries
    print(f"[*] Loading paddock boundary shapefile: {args.shapefile}")
    gdf = gpd.read_file(args.shapefile)
    gdf_wgs84 = gdf.to_crs("EPSG:4326")
    geom_multi = gdf_wgs84.geometry.unary_union
    bbox_coords = list(geom_multi.bounds)
    
    # 2. Load Growth Spike and Frequency Rasters
    print(f"[*] Loading targeted rain-pulse growth spike: {args.max_spike}")
    spike_raw, transform_ref, crs_ref = read_raster_raw(
        args.max_spike, cmap_name='YlGn', vmin=0.0, vmax=0.40
    )
    
    print(f"[*] Loading targeted rain-pulse growth frequency: {args.growth_freq}")
    freq_raw, _, _ = read_raster_raw(
        args.growth_freq, cmap_name='plasma', vmin=0.0, vmax=1.0
    )
    
    da_spike_template = xr.DataArray(
        spike_raw,
        coords={'y': np.arange(spike_raw.shape[0]), 'x': np.arange(spike_raw.shape[1])},
        dims=['y', 'x']
    )
    # Open spike raster via rioxarray for geographic template matching
    da_spike_rio = rioxarray.open_rasterio(args.max_spike)
    if da_spike_rio.shape[0] > 1:
        # Convert multi-band template to single-band coordinates
        da_spike_template = xr.DataArray(
            spike_raw,
            coords={'y': da_spike_rio.y, 'x': da_spike_rio.x},
            dims=['y', 'x']
        ).rio.write_crs(crs_ref)
    else:
        da_spike_template = da_spike_rio.squeeze(axis=0).rio.write_crs(crs_ref)

    # 3. Query AWS Earth Search STAC API
    print("[*] Connecting to Earth Search STAC API...")
    catalog = pystac_client.Client.open("https://earth-search.aws.element84.com/v1")
    
    search_range = f"{args.start_year}-08-01/{args.end_year}-11-15"
    print(f"    [*] Querying Sentinel-2 scenes covering {search_range} with cloud cover < 10%...")
    
    search = catalog.search(
        collections=["sentinel-2-l2a"],
        bbox=bbox_coords,
        datetime=search_range,
        query={"eo:cloud_cover": {"lt": 10}}
    )
    
    max_retries = 3
    items = None
    for attempt in range(max_retries):
        try:
            items = search.item_collection()
            break
        except Exception as e:
            if attempt == max_retries - 1:
                raise e
            print(f"    [!] STAC query failed ({e}). Retrying in 5 seconds...")
            time.sleep(5)
            
    # Filter to dry-season months (August, September, October, November)
    dry_season_items = []
    for it in items:
        dt_str = it.properties["datetime"]
        month = int(dt_str[5:7])
        if month in [8, 9, 10, 11]:
            dry_season_items.append(it)
            
    print(f"[+] Found {len(dry_season_items)} clean dry-season scenes across {args.start_year}-{args.end_year}.")
    
    # 4. Filter Baseline and Curing Scenes
    sorted_baseline_items = sorted(dry_season_items, key=lambda x: x.properties["eo:cloud_cover"])[:6]
    
    # Select earliest and latest scenes of the most recent year in the query for curing rate
    most_recent_year = str(args.end_year)
    items_recent = [it for it in dry_season_items if most_recent_year in it.properties["datetime"]]
    items_recent_sorted = sorted(items_recent, key=lambda x: x.properties["datetime"])
    
    if len(items_recent_sorted) >= 2:
        earliest_recent = items_recent_sorted[0]
        latest_recent = items_recent_sorted[-1]
        curing_candidates = [earliest_recent, latest_recent]
    else:
        earliest_recent = None
        latest_recent = None
        curing_candidates = []
        print(f"    [!] Warning: Only found {len(items_recent_sorted)} scenes in {most_recent_year}. Curing rate will be skipped.")
        
    unique_items = []
    seen_ids = set()
    for item in sorted_baseline_items + curing_candidates:
        if item.id not in seen_ids:
            unique_items.append(item)
            seen_ids.add(item.id)
            
    msavi_dict = {}
    ndre_dict = {}
    ndwi_dict = {}
    
    print(f"[*] Processing {len(unique_items)} unique Sentinel-2 scenes end-to-end...")
    for idx, item in enumerate(unique_items):
        try:
            date_str = item.properties["datetime"][:10]
            cloud_cover = item.properties["eo:cloud_cover"]
            print(f"    [*] Loading scene {idx+1}/{len(unique_items)} ({date_str}, cloud cover: {cloud_cover:.1f}%)...", flush=True)
            
            # Open assets
            red_url = item.assets["red"].href
            nir_url = item.assets["nir"].href
            scl_url = item.assets["scl"].href
            b05_url = item.assets["rededge1"].href
            b11_url = item.assets["swir16"].href
            
            da_red = rioxarray.open_rasterio(red_url)
            da_nir = rioxarray.open_rasterio(nir_url)
            da_scl = rioxarray.open_rasterio(scl_url)
            da_b05 = rioxarray.open_rasterio(b05_url)
            da_b11 = rioxarray.open_rasterio(b11_url)
            
            # Reproject boundaries to native scene CRS
            s2_crs = da_red.rio.crs
            gdf_s2 = gdf.to_crs(s2_crs)
            s2_bbox = gdf_s2.total_bounds
            geom_s2 = gdf_s2.geometry.unary_union.buffer(args.buffer_pixels * 10.0)
            
            # Crop via HTTP range requests
            da_red_cropped = da_red.rio.clip_box(*s2_bbox).rio.clip([geom_s2], crs=s2_crs, all_touched=True)
            da_nir_cropped = da_nir.rio.clip_box(*s2_bbox).rio.clip([geom_s2], crs=s2_crs, all_touched=True)
            da_scl_cropped = da_scl.rio.clip_box(*s2_bbox).rio.clip([geom_s2], crs=s2_crs, all_touched=True)
            da_b05_cropped = da_b05.rio.clip_box(*s2_bbox).rio.clip([geom_s2], crs=s2_crs, all_touched=True)
            da_b11_cropped = da_b11.rio.clip_box(*s2_bbox).rio.clip([geom_s2], crs=s2_crs, all_touched=True)
            
            # Convert reflectances (divide by 10000)
            red_vals = da_red_cropped.values.astype(float).squeeze(axis=0) / 10000.0
            nir_vals = da_nir_cropped.values.astype(float).squeeze(axis=0) / 10000.0
            
            # Resample 20m bands to 10m grid
            da_scl_resampled = da_scl_cropped.rio.reproject_match(da_red_cropped, resampling=Resampling.nearest)
            da_b05_resampled = da_b05_cropped.rio.reproject_match(da_red_cropped, resampling=Resampling.bilinear)
            da_b11_resampled = da_b11_cropped.rio.reproject_match(da_red_cropped, resampling=Resampling.bilinear)
            
            scl_vals = da_scl_resampled.values.squeeze(axis=0)
            b05_vals = da_b05_resampled.values.astype(float).squeeze(axis=0) / 10000.0
            b11_vals = da_b11_resampled.values.astype(float).squeeze(axis=0) / 10000.0
            
            # Mask clouds and sensor artifacts
            cloud_mask = np.isin(scl_vals, [0, 1, 3, 8, 9, 10])
            red_vals[cloud_mask] = np.nan
            nir_vals[cloud_mask] = np.nan
            b05_vals[cloud_mask] = np.nan
            b11_vals[cloud_mask] = np.nan
            
            # Calculate indices
            sqrt_term = (2.0 * nir_vals + 1.0)**2 - 8.0 * (nir_vals - red_vals)
            sqrt_term = np.clip(sqrt_term, 0, None)
            msavi = (2.0 * nir_vals + 1.0 - np.sqrt(sqrt_term)) / 2.0
            
            ndre = (nir_vals - b05_vals) / np.clip(nir_vals + b05_vals, 1e-5, None)
            ndwi = (nir_vals - b11_vals) / np.clip(nir_vals + b11_vals, 1e-5, None)
            
            # Wrap as spatial rasters
            da_msavi = xr.DataArray(msavi, coords={'y': da_red_cropped.y, 'x': da_red_cropped.x}, dims=['y', 'x']).rio.write_crs(s2_crs)
            da_ndre = xr.DataArray(ndre, coords={'y': da_red_cropped.y, 'x': da_red_cropped.x}, dims=['y', 'x']).rio.write_crs(s2_crs)
            da_ndwi = xr.DataArray(ndwi, coords={'y': da_red_cropped.y, 'x': da_red_cropped.x}, dims=['y', 'x']).rio.write_crs(s2_crs)
            
            # Align immediately to template matching project grid
            msavi_dict[item.id] = da_msavi.rio.reproject_match(da_spike_template)
            ndre_dict[item.id] = da_ndre.rio.reproject_match(da_spike_template)
            ndwi_dict[item.id] = da_ndwi.rio.reproject_match(da_spike_template)
            
            print(f"    [+] Loaded scene {idx+1}/{len(unique_items)}", flush=True)
            
        except Exception as e:
            print(f"    [!] Error loading scene {item.id}: {e}", flush=True)
            continue
            
    # Calculate statistics
    loaded_baseline = [it for it in sorted_baseline_items if it.id in msavi_dict]
    if len(loaded_baseline) < 2:
        raise RuntimeError("Too few successful scenes loaded to build a reliable baseline.")
        
    print("[*] Computing climatological baseline statistics...")
    msavi_stack = xr.concat([msavi_dict[it.id] for it in loaded_baseline], dim='time')
    mean_msavi = msavi_stack.mean(dim='time', skipna=True).values.squeeze()
    min_msavi = msavi_stack.min(dim='time', skipna=True).values.squeeze()
    std_msavi = msavi_stack.std(dim='time', skipna=True).values.squeeze()
    
    ndre_stack = xr.concat([ndre_dict[it.id] for it in loaded_baseline], dim='time')
    mean_ndre = ndre_stack.mean(dim='time', skipna=True).values.squeeze()
    
    ndwi_stack = xr.concat([ndwi_dict[it.id] for it in loaded_baseline], dim='time')
    mean_ndwi = ndwi_stack.mean(dim='time', skipna=True).values.squeeze()
    
    # Dry-season senescence verification: confirms grass cures down during dry spells
    print("[*] Computing dry-season senescence mask (min dry-season MSAVI < 0.22)...")
    curing_mask = np.isnan(min_msavi) | (min_msavi < 0.22)
        
    # Apply 3-Tier Classification Decision Tree
    print("[*] Applying Decision Tree classification rules...")
    tiers = np.zeros_like(spike_raw, dtype=np.uint8)
    
    valid_mask = ~np.isnan(spike_raw)
    valid_msavi = ~np.isnan(mean_msavi) & ~np.isnan(std_msavi) & ~np.isnan(mean_ndre) & ~np.isnan(mean_ndwi)
    
    # Tier 3: Native Refugia (Evergreen Tree Canopy)
    tier3_mask = (
        valid_mask & valid_msavi &
        (mean_msavi >= 0.20) &
        (mean_ndre >= 0.15) &
        (mean_ndwi >= -0.22) &
        (std_msavi < 0.06) &
        (spike_raw < 0.12)
    )
    
    # Tier 2: Heavy Infestation (Core containment)
    tier2_mask = valid_mask & (spike_raw >= 0.20) & (freq_raw >= 0.30) & (~tier3_mask) & curing_mask
    
    # Tier 1: Light Infestation (Active eradication)
    tier1_mask = valid_mask & ((spike_raw >= 0.12) & (freq_raw >= 0.10)) & (~tier3_mask) & (~tier2_mask) & curing_mask
    
    tiers[tier1_mask] = 1
    tiers[tier2_mask] = 2
    tiers[tier3_mask] = 3
    
    # 5. Save Classified Outputs
    tif_out = os.path.join(args.output_dir, "buffel_management_tiers.tif")
    png_out = os.path.join(args.output_dir, "buffel_management_tiers.png")
    html_out = os.path.join(args.output_dir, "keewong_management_tiers_map.html")
    csv_out = os.path.join(args.output_dir, "management_tier_statistics.csv")
    
    # Save GeoTIFF
    da_tiers = xr.DataArray(
        tiers,
        coords={'y': da_spike_template.y, 'x': da_spike_template.x},
        dims=['y', 'x']
    ).rio.write_crs(crs_ref)
    
    # Clean transparent out-of-boundary pixels (set NaN background to 0/transparent)
    da_tiers_masked = da_tiers.where(valid_mask, 0)
    da_tiers_masked.rio.to_raster(tif_out, dtype=np.uint8)
    print(f"[+] Saved classified GeoTIFF to: {tif_out}")
    
    # Save colorized PNG (reprojected to WGS84 for Leaflet layer overlay)
    da_tiers_wgs84 = da_tiers_masked.rio.reproject("EPSG:4326", resampling=Resampling.nearest)
    wgs84_data = da_tiers_wgs84.values.squeeze()
    
    h, w = wgs84_data.shape
    png_rgba = np.zeros((h, w, 4), dtype=np.uint8)
    
    # Map colors: 1->Green, 2->Red, 3->Cyan, 0->transparent
    png_rgba[wgs84_data == 1] = [44, 160, 44, 215]   # Green (Eradication)
    png_rgba[wgs84_data == 2] = [214, 39, 40, 215]   # Red (Containment)
    png_rgba[wgs84_data == 3] = [0, 255, 255, 215]   # Cyan (Refugia)
    
    from PIL import Image
    im = Image.fromarray(png_rgba)
    im.save(png_out)
    print(f"[+] Saved colorized PNG overlay to: {png_out}")
    
    # Save Leaflet HTML dashboard
    wgs84_bounds = list(da_tiers_wgs84.rio.bounds())
    generate_leaflet_map(da_tiers_wgs84, gdf_wgs84, html_out, wgs84_bounds)
    print(f"[+] Saved Leaflet HTML map dashboard to: {html_out}")
    
    # Compute Statistics
    pixel_area_ha = (abs(transform_ref[0] * transform_ref[4])) / 10000.0
    
    t1_area = np.sum(tiers == 1) * pixel_area_ha
    t2_area = np.sum(tiers == 2) * pixel_area_ha
    t3_area = np.sum(tiers == 3) * pixel_area_ha
    bg_area = np.sum((tiers == 0) & valid_mask) * pixel_area_ha
    total_area = np.sum(valid_mask) * pixel_area_ha
    
    stats = pd.DataFrame([
        {"Tier": "Tier 1 (Eradication)", "Area_ha": round(t1_area, 2), "Pct": round((t1_area/total_area)*100.0, 1)},
        {"Tier": "Tier 2 (Containment)", "Area_ha": round(t2_area, 2), "Pct": round((t2_area/total_area)*100.0, 1)},
        {"Tier": "Tier 3 (Refugia)", "Area_ha": round(t3_area, 2), "Pct": round((t3_area/total_area)*100.0, 1)},
        {"Tier": "Background/Soil", "Area_ha": round(bg_area, 2), "Pct": round((bg_area/total_area)*100.0, 1)},
        {"Tier": "Total Area", "Area_ha": round(total_area, 2), "Pct": 100.0}
    ])
    stats.to_csv(csv_out, index=False)
    print(f"[+] Saved CSV statistics to: {csv_out}")
    
    print("\n==========================================")
    print("         MANAGEMENT TIER AREAS")
    print("==========================================")
    print(f"Total Area:           {total_area:.2f} ha")
    print(f"Tier 1 (Eradication): {t1_area:.2f} ha ({t1_area/total_area*100.0:.1f}%)")
    print(f"Tier 2 (Containment): {t2_area:.2f} ha ({t2_area/total_area*100.0:.1f}%)")
    print(f"Tier 3 (Refugia):     {t3_area:.2f} ha ({t3_area/total_area*100.0:.1f}%)")
    print("==========================================\n")
    print("[+] 3-TIER LAYER GENERATION COMPLETE!")
    
if __name__ == "__main__":
    main()
