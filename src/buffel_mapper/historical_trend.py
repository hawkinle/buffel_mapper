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
        geom = gdf_wgs84.geometry.iloc[0]
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

def save_colorized_tiff(data_2d, vmin, vmax, cmap_name, out_path, crs, transform):
    """
    Save a 2D single-band float array as a colorized 3-band RGB GeoTIFF.
    """
    cmap = plt.get_cmap(cmap_name)
    norm = np.clip((data_2d - vmin) / (vmax - vmin), 0, 1)
    rgba = cmap(norm)
    rgb = rgba[..., :3]
    rgb_band_first = np.transpose(rgb, (2, 0, 1)).astype(np.float32)
    
    meta = {
        'driver': 'GTiff',
        'dtype': 'float32',
        'nodata': None,
        'width': rgb.shape[1],
        'height': rgb.shape[0],
        'count': 3,
        'crs': crs,
        'transform': transform
    }
    with rasterio.open(out_path, 'w', **meta) as dst:
        dst.write(rgb_band_first)

def run_historical_trend(kml_path, start_year, end_year, max_cloud, out_dir, threshold=0.12, buffer_pixels=0):
    """
    Compute annual NDVI standard deviation maps for each year and compile trend statistics.
    """
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
    
    print(f"[*] Output directory set to: {out_dir}")
    print(f"[*] Bounding Box: {bbox_coords} (WGS84)")
    print(f"[*] Year range: {start_year} to {end_year}")
    
    # 1. Connect to Planetary Computer STAC catalog
    print("[*] Connecting to Planetary Computer STAC API...")
    catalog = pystac_client.Client.open(
        "https://planetarycomputer.microsoft.com/api/stac/v1",
        modifier=planetary_computer.sign_inplace
    )
    
    geom = box(*bbox_coords)
    gdf = gpd.GeoDataFrame(geometry=[geom], crs="EPSG:4326")
    
    results = []
    
    leaflet_bounds = None
    crs_ref = None
    transform_ref = None
    
    # Process year-by-year
    for year in range(start_year, end_year + 1):
        print(f"\n==================== PROCESSING YEAR {year} ====================")
        search = catalog.search(
            collections=["sentinel-2-l2a"],
            bbox=bbox_coords,
            datetime=f"{year}-01-01/{year}-12-31",
            query={"eo:cloud_cover": {"lt": max_cloud}}
        )
        
        items = search.item_collection()
        if len(items) == 0:
            print(f"[!] No scenes found for year {year}. Skipping.")
            continue
            
        print(f"[+] Found {len(items)} scenes for year {year}. Loading NDVI...")
        
        ndvi_list = []
        
        # Sort scenes chronologically
        sorted_items = sorted(items, key=lambda x: x.properties["datetime"])
        
        for i, item in enumerate(sorted_items):
            try:
                red_url = item.assets["B04"].href
                nir_url = item.assets["B08"].href
                scl_url = item.assets["SCL"].href
                
                # Load and crop Red
                da_red = rioxarray.open_rasterio(red_url)
                gdf_proj = gdf.to_crs(da_red.rio.crs)
                proj_bbox = gdf_proj.total_bounds
                clip_geom, clip_crs = get_clip_geometry(polygon, da_red.rio.crs, buffer_pixels)
                red_crop = da_red.rio.clip_box(*proj_bbox).rio.clip([clip_geom], crs=clip_crs, all_touched=True)
                red_val = red_crop.where(red_crop != da_red.rio.nodata).astype(np.float32) / 10000.0
                
                # Load and crop NIR
                da_nir = rioxarray.open_rasterio(nir_url)
                nir_crop = da_nir.rio.clip_box(*proj_bbox).rio.clip([clip_geom], crs=clip_crs, all_touched=True)
                nir_val = nir_crop.where(nir_crop != da_nir.rio.nodata).astype(np.float32) / 10000.0
                
                # Load and crop SCL
                da_scl = rioxarray.open_rasterio(scl_url)
                scl_crop = da_scl.rio.clip_box(*proj_bbox).rio.clip([clip_geom], crs=clip_crs, all_touched=True)
                scl_resampled = scl_crop.rio.reproject_match(red_val)
                scl_val = scl_resampled.squeeze()
                
                # Mask out cloud/shadow
                invalid_scl = [0, 1, 3, 8, 9, 10, 11]
                valid_mask = ~np.isin(scl_val.values, invalid_scl)
                
                red_masked = red_val.where(valid_mask)
                nir_masked = nir_val.where(valid_mask)
                
                ndvi = (nir_masked - red_masked) / (nir_masked + red_masked)
                ndvi_list.append(ndvi)
                
                # Print progress bar every 5 scenes or for the last scene
                if (i + 1) % 5 == 0 or (i + 1) == len(sorted_items):
                    percent_done = int((i + 1) / len(sorted_items) * 100)
                    bar_len = 20
                    filled_len = int(bar_len * (i + 1) // len(sorted_items))
                    bar = '=' * filled_len + '-' * (bar_len - filled_len)
                    print(f"    [{bar}] {percent_done}% - Loaded scene {i+1}/{len(sorted_items)}", flush=True)
                
                # Capture metadata references on the first successful load
                if crs_ref is None:
                    crs_ref = ndvi.rio.crs
                    transform_ref = ndvi.rio.transform()
                    da_gps = ndvi.rio.reproject("EPSG:4326")
                    bounds = da_gps.rio.bounds()
                    leaflet_bounds = [[bounds[1], bounds[0]], [bounds[3], bounds[2]]]
                    
            except Exception as e:
                # Skip individual scene failures silently
                continue
                
        if len(ndvi_list) < 5:
            print(f"[!] Too few successful scenes for year {year} ({len(ndvi_list)}). Skipping statistics.")
            continue
            
        print(f"[*] Calculating standard deviation for year {year} ({len(ndvi_list)} scenes)...")
        ndvi_stack = xr.concat(ndvi_list, dim='time')
        ndvi_std = ndvi_stack.std(dim='time', skipna=True)
        ndvi_std_2d = ndvi_std.values.squeeze()
        
        # Save colorized probability GeoTIFF for this year
        year_tif_path = os.path.join(out_dir, f"buffel_probability_{year}.tif")
        save_colorized_tiff(ndvi_std_2d, 0.05, 0.20, 'viridis', year_tif_path, crs_ref, transform_ref)
        
        # Save visual PNG preview for Leaflet Overlay
        year_png_path = os.path.join(out_dir, f"buffel_probability_{year}.png")
        plt.imsave(year_png_path, ndvi_std_2d, cmap="viridis", vmin=0.05, vmax=0.20)
        
        # Calculate area statistics (1 pixel = 0.01 ha)
        pixel_area_ha = 0.01
        valid_mask = ~np.isnan(ndvi_std_2d)
        total_pixels = np.sum(valid_mask)
        total_area_ha = total_pixels * pixel_area_ha
        
        infested_mask = (ndvi_std_2d >= threshold) & valid_mask
        infested_pixels = np.sum(infested_mask)
        infested_area_ha = infested_pixels * pixel_area_ha
        
        infestation_rate = (infested_area_ha / total_area_ha) * 100 if total_area_ha > 0 else 0
        
        print(f"[+] Year {year} Results: Total Area = {total_area_ha:.2f} ha, Infested = {infested_area_ha:.2f} ha ({infestation_rate:.1f}%)")
        
        results.append({
            'Year': year,
            'Total_Area_Ha': total_area_ha,
            'Infested_Area_Ha': infested_area_ha,
            'Infestation_Rate_Pct': infestation_rate,
            'Image_File': f"buffel_probability_{year}.png"
        })
        
    if not results:
        print("[!] No years were successfully processed.")
        return
        
    # 5. Save results to CSV
    df = pd.DataFrame(results)
    csv_path = os.path.join(out_dir, "historical_infestation_trend.csv")
    df.to_csv(csv_path, index=False)
    print(f"\n[+] Saved historical trend data to: {csv_path}")
    
    # 6. Generate Historical Trend Line Plot
    plt.figure(figsize=(11, 6))
    plt.plot(df['Year'], df['Infested_Area_Ha'], marker='o', color='#ef4444', linewidth=2.5, markersize=8)
    plt.title('Historical Buffel Grass Infestation Trend (Keewong)', fontsize=14, fontweight='bold', pad=15)
    plt.xlabel('Year', fontsize=12)
    plt.ylabel('Estimated Infestation (Hectares)', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xticks(df['Year'])
    
    # Add labels on points
    for idx, row in df.iterrows():
        plt.text(row['Year'], row['Infested_Area_Ha'] + (max(df['Infested_Area_Ha'])*0.02), 
                 f"{row['Infested_Area_Ha']:.1f} ha", ha='center', fontsize=9, fontweight='bold')
                 
    plot_path = os.path.join(out_dir, "historical_trend_chart.png")
    plt.tight_layout()
    plt.savefig(plot_path, dpi=150)
    plt.close()
    print(f"[+] Saved historical trend line plot: {plot_path}")
    
    # 7. Generate Leaflet HTML Map with Year Selector
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
    
    # Generate Javascript layer dictionary
    js_layers_init = ""
    js_layers_dict = ""
    for r in results:
        y = r['Year']
        img = r['Image_File']
        js_layers_init += f"    var prob_{y} = L.imageOverlay('{img}', {leaflet_bounds}, {{ opacity: 0.65, interactive: false }});\n"
        # Make the most recent year default visible
        if y == results[-1]['Year']:
            js_layers_init += f"    prob_{y}.addTo(map);\n"
        js_layers_dict += f"        'Infestation {y}': prob_{y},\n"
        
    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <title>Keewong Property - Historical Weed Map (2016-2025)</title>
    <meta charset="utf-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <!-- Leaflet CSS -->
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
        /* Info Panel */
        .info-panel {{
            padding: 15px;
            background: rgba(30, 30, 36, 0.95);
            box-shadow: 0 4px 15px rgba(0,0,0,0.5);
            border-radius: 8px;
            border: 1px solid #444;
            color: #f3f3f3;
            max-width: 320px;
        }}
        .info-panel h2 {{
            margin: 0 0 10px 0;
            font-size: 16px;
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
        /* Legend */
        .legend {{
            padding: 10px;
            background: rgba(30, 30, 36, 0.95);
            border: 1px solid #444;
            border-radius: 8px;
            line-height: 18px;
            color: #f3f3f3;
            font-size: 12px;
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

<!-- Leaflet JS -->
<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
<script>
    // Initialize map
    var map = L.map('map').setView([{center_lat}, {center_lon}], 14);

    // Add Esri World Imagery Basemap
    var satellite = L.tileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{{z}}/{{y}}/{{x}}', {{
        attribution: 'Tiles &copy; Esri &mdash; Source: Esri'
    }}).addTo(map);

    // Property Boundary Polygon (from KML)
    var boundaryCoords = {leaflet_coords};
    var boundary = L.polygon(boundaryCoords, {{
        color: '#38bdf8',
        fillColor: 'transparent',
        weight: 3,
        opacity: 0.9
    }}).addTo(map);

    // Load all annual overlay layers
{js_layers_init}

    // Define Layer Controls
    var basemaps = {{
        "Satellite Imagery": satellite
    }};
    
    var overlayMaps = {{
        "Property Boundary": boundary,
{js_layers_dict}
    }};
    
    L.control.layers(basemaps, overlayMaps, {{ collapsed: false, position: 'topright' }}).addTo(map);

    // Add Area Statistics Control Panel
    var info = L.control({{ position: 'bottomleft' }});
    info.onAdd = function (map) {{
        var div = L.DomUtil.create('div', 'info-panel');
        var html = '<h2>Keewong Property - Historical Stats</h2>' +
            '<div class="stat-row"><span>Total Property:</span><span class="stat-val">{total_area_ha:.1f} ha</span></div>' +
            '<div style="font-size: 12px; font-weight: bold; border-top: 1px solid #444; margin-top: 8px; padding-top: 8px; color: #4ade80;">Infestation Area by Year:</div>';
        var yearData = {" + ".join([f"\'<div class=\\\"stat-row\\\"><span>Year {int(r['Year'])}:</span><span class=\\\"stat-highlight\\\">{r['Infested_Area_Ha']:.1f} ha ({r['Infestation_Rate_Pct']:.1f}%)</span></div>\'" for r in results])};
            
        div.innerHTML = html + yearData +
            '<div class="stat-row" style="font-size: 11px; color:#aaa; margin-top:10px; border-top:1px solid #444; padding-top:5px;">' +
            '<span>Analysis Method:</span><span>SCL Cloud Masking + NDVI SD</span></div>';
        return div;
    }};
    info.addTo(map);

    // Add Legend Control Panel
    var legend = L.control({{ position: 'bottomright' }});
    legend.onAdd = function (map) {{
        var div = L.DomUtil.create('div', 'legend');
        div.innerHTML = '<h3 style="margin:0 0 5px 0; font-size:12px; color:#4ade80;">Weed Volatility</h3>' +
            '<i style="background: #fde725"></i> High Volatility (Active Growth)<br>' +
            '<i style="background: #21918c"></i> Moderate Growth<br>' +
            '<i style="background: #440154"></i> Stable / Bare Soil<br>';
        return div;
    }};
    legend.addTo(map);

    // Adjust map viewport
    map.fitBounds(boundary.getBounds());
</script>

</body>
</html>
"""
    
    html_path = os.path.join(out_dir, "keewong_historical_weed_map.html")
    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    print(f"[+] Saved historical Leaflet HTML map: {html_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Decade Historical Weed Trend Analyzer")
    parser.add_argument("--kml", type=str, required=True, help="Path to KML file")
    parser.add_argument("--start-year", type=int, default=2016, help="Start Year (2016-)")
    parser.add_argument("--end-year", type=int, default=2025, help="End Year")
    parser.add_argument("--max-cloud", type=float, default=10.0, help="Max cloud cover %")
    parser.add_argument("--out-dir", type=str, default="output_historical", help="Output directory")
    parser.add_argument("--threshold", type=float, default=0.12, help="NDVI Std Dev threshold")
    parser.add_argument("--buffer-pixels", type=int, default=0, help="Buffer boundaries by N pixels (1 pixel = 10m)")
    
    args = parser.parse_args()
    run_historical_trend(args.kml, args.start_year, args.end_year, args.max_cloud, args.out_dir, args.threshold, args.buffer_pixels)
