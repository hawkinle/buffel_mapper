#!/usr/bin/env python3
import os
import argparse
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import rioxarray
import shapely.geometry
import json

def get_clip_geometry(polygon, raster_crs, buffer_pixels=2):
    gdf_poly = gpd.GeoDataFrame(geometry=[polygon], crs="EPSG:4326")
    gdf_poly_proj = gdf_poly.to_crs(raster_crs)
    if buffer_pixels > 0:
        buffer_meters = buffer_pixels * 10.0
        buffered_geom = gdf_poly_proj.geometry.buffer(buffer_meters).iloc[0]
        return buffered_geom
    else:
        return gdf_poly_proj.geometry.iloc[0]

def rgb_to_float(rgb_array, cmap_name, vmin, vmax):
    """
    Reconstruct raw float data from a colorized 3-band RGB array using pure NumPy.
    """
    # rgb_array shape: (3, h, w)
    h, w = rgb_array.shape[1], rgb_array.shape[2]
    rgb_3d = np.transpose(rgb_array, (1, 2, 0))
    flat_rgb = rgb_3d.reshape(-1, 3)
    
    # Identify background/NaN: where RGB is all zero, or has NaNs
    is_bg = np.all(flat_rgb < 0.01, axis=1) | np.any(np.isnan(flat_rgb), axis=1)
    
    # Generate colormap lookup table
    x = np.linspace(0, 1, 1000)
    cmap = plt.get_cmap(cmap_name)
    lut = cmap(x)[..., :3] # shape (1000, 3)
    
    raw_vals = np.full(flat_rgb.shape[0], np.nan)
    non_bg_rgb = flat_rgb[~is_bg]
    
    if len(non_bg_rgb) > 0:
        # Loop over 1000 colors in LUT to find the nearest match for each pixel
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

def clip_and_reconstruct(tif_path, buffered_geom, raster_crs, cmap_name, vmin, vmax):
    da = rioxarray.open_rasterio(tif_path)
    clipped = da.rio.clip([buffered_geom], crs=raster_crs, all_touched=True)
    vals = clipped.values
    if vals.shape[0] == 1:
        raw = vals[0].astype(float)
        nodata = da.rio.nodata
        if nodata is not None and not np.isnan(nodata):
            raw[raw == nodata] = np.nan
        return raw
    return rgb_to_float(vals, cmap_name, vmin, vmax)

def assign_dynamic_paddock_name(lat, lon, bounds, idx, used_names):
    """Assign a descriptive compass name based on relative position within property bounds."""
    min_lon, min_lat, max_lon, max_lat = bounds
    lat_span = max(max_lat - min_lat, 1e-6)
    lon_span = max(max_lon - min_lon, 1e-6)
    rel_y = (lat - min_lat) / lat_span  # 0 = South, 1 = North
    rel_x = (lon - min_lon) / lon_span  # 0 = West, 1 = East

    if rel_y > 0.55 and rel_x > 0.55:
        base = "Central-East Paddock"
    elif rel_y < 0.25 and rel_x > 0.70:
        base = "South-East Paddock"
    elif rel_y < 0.45 and rel_x > 0.55:
        base = "South-Central Paddock"
    elif rel_y < 0.55 and rel_x < 0.45:
        base = "Western Paddock"
    elif rel_y >= 0.55 and rel_x <= 0.55:
        base = "North-Western Paddock"
    else:
        ns = "North" if rel_y >= 0.6 else ("South" if rel_y <= 0.4 else "Central")
        ew = "East" if rel_x >= 0.6 else ("West" if rel_x <= 0.4 else "Central")
        base = f"{ns}-{ew} Paddock" if ns != ew else f"{ns} Paddock"

    name = base
    if name in used_names:
        name = f"{base} ({idx + 1})"
    used_names.add(name)
    return name

def main():
    parser = argparse.ArgumentParser(description="Paddock-by-Paddock Buffel Infestation Breakdown")
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
        
    out_dir = os.path.join(project_dir, "paddock_breakdown")
    os.makedirs(out_dir, exist_ok=True)
    
    print("=" * 60)
    print(f"[*] RUNNING PADDOCK-BY-PADDOCK BREAKDOWN FOR PROJECT: {args.project}")
    print(f"[*] Shapefile: {shapefile_path}")
    print(f"[*] Output Directory: {out_dir}")
    print("=" * 60)
    
    # 1. Load shapefile and extract individual paddocks
    gdf = gpd.read_file(shapefile_path)
    gdf_wgs84 = gdf.to_crs("EPSG:4326")
    geom_multi = gdf_wgs84.geometry.unary_union
    prop_bounds = geom_multi.bounds
    
    if not isinstance(geom_multi, shapely.geometry.MultiPolygon):
        print("[!] Warning: Geometry is not a MultiPolygon. Treating as single paddock.")
        paddock_polygons = [geom_multi]
        paddock_polygons_proj = [gdf.geometry.unary_union]
    else:
        paddock_polygons = list(geom_multi.geoms)
        paddock_polygons_proj = list(gdf.geometry.unary_union.geoms)
        
    paddocks = []
    used_names = set()
    for idx, (poly_wgs84, poly_proj) in enumerate(zip(paddock_polygons, paddock_polygons_proj)):
        # Calculate nominal area (pre-buffer)
        nominal_area = poly_proj.area / 10000.0
        lat, lon = poly_wgs84.centroid.y, poly_wgs84.centroid.x
        name = assign_dynamic_paddock_name(lat, lon, prop_bounds, idx, used_names)
            
        paddocks.append({
            'index': idx,
            'name': name,
            'poly_wgs84': poly_wgs84,
            'poly_proj': poly_proj,
            'nominal_area_ha': nominal_area,
            'centroid': (lat, lon)
        })
        
    # Sort paddocks from largest to smallest area
    paddocks = sorted(paddocks, key=lambda x: x['nominal_area_ha'], reverse=True)
    
    # Load targeted rasters reference for CRS
    max_spike_path = os.path.join(project_dir, "targeted", "buffel_max_growth_spike.tif")
    freq_path = os.path.join(project_dir, "targeted", "buffel_growth_frequency.tif")
    tiers_path = os.path.join(project_dir, "tiers_layer", "buffel_management_tiers.tif")
    
    if not os.path.exists(max_spike_path) or not os.path.exists(freq_path):
        raise FileNotFoundError("Targeted analysis rasters (.tif) not found in targeted/ folder.")
        
    da_ref = rioxarray.open_rasterio(max_spike_path)
    raster_crs = da_ref.rio.crs
    
    summary_data = []
    historical_trends = {}
    
    # Years to scan
    years = [y for y in range(2016, 2026) if y != 2022] # skip 2022
    
    for pad in paddocks:
        print(f"\n[*] Processing Paddock: {pad['name']} ({pad['nominal_area_ha']:.2f} ha)...")
        buffered_geom = get_clip_geometry(pad['poly_wgs84'], raster_crs, args.buffer_pixels)
        
        # 1. Clip and load targeted values
        spike_vals = clip_and_reconstruct(max_spike_path, buffered_geom, raster_crs, 'YlGn', 0.0, 0.40)
        freq_vals = clip_and_reconstruct(freq_path, buffered_geom, raster_crs, 'plasma', 0.0, 1.0)
        
        valid_mask = ~np.isnan(spike_vals)
        total_buffered_ha = np.sum(valid_mask) * 0.01
        
        # High confidence core condition (calibrated for 1-band P90 MSAVI vs legacy 3-band Max NDVI)
        is_p90_msavi = (da_ref.shape[0] == 1)
        s_thresh = 0.105 if is_p90_msavi else 0.20
        f_thresh = 0.050 if is_p90_msavi else 0.30
        high_conf_mask = (spike_vals >= s_thresh) & (freq_vals >= f_thresh) & valid_mask
        infested_ha = np.sum(high_conf_mask) * 0.01
        infestation_rate = (infested_ha / total_buffered_ha * 100) if total_buffered_ha > 0 else 0.0
        
        # Load 3-Tier breakdown if tiers_layer exists
        t1_ha, t2_ha, t3_ha = np.nan, np.nan, np.nan
        t1_pct, t2_pct, t3_pct = np.nan, np.nan, np.nan
        if os.path.exists(tiers_path):
            tier_vals = clip_and_reconstruct(tiers_path, buffered_geom, raster_crs, 'viridis', 0, 3)
            t1_ha = np.sum((tier_vals == 1) & valid_mask) * 0.01
            t2_ha = np.sum((tier_vals == 2) & valid_mask) * 0.01
            t3_ha = np.sum((tier_vals == 3) & valid_mask) * 0.01
            if total_buffered_ha > 0:
                t1_pct = (t1_ha / total_buffered_ha) * 100.0
                t2_pct = (t2_ha / total_buffered_ha) * 100.0
                t3_pct = (t3_ha / total_buffered_ha) * 100.0
        
        print(f"    [+] Buffered Area: {total_buffered_ha:.2f} ha")
        print(f"    [+] Targeted Core Infested: {infested_ha:.2f} ha ({infestation_rate:.1f}%)")
        if not np.isnan(t1_ha):
            print(f"    [+] Tiers -> T1 (Eradicate): {t1_ha:.2f} ha ({t1_pct:.1f}%) | T2 (Contain): {t2_ha:.2f} ha ({t2_pct:.1f}%) | T3 (Refugia): {t3_ha:.2f} ha ({t3_pct:.1f}%)")
        
        # 2. Historical stats
        pad_trend = {}
        for year in years:
            year_tif_path = os.path.join(project_dir, "historical", f"buffel_probability_{year}.tif")
            if os.path.exists(year_tif_path):
                std_vals = clip_and_reconstruct(year_tif_path, buffered_geom, raster_crs, 'viridis', 0.05, 0.20)
                valid_std = ~np.isnan(std_vals)
                total_std_ha = np.sum(valid_std) * 0.01
                
                # Volatility condition: std >= 0.12
                infested_std_mask = (std_vals >= 0.12) & valid_std
                infested_std_ha = np.sum(infested_std_mask) * 0.01
                infested_std_pct = (infested_std_ha / total_std_ha * 100) if total_std_ha > 0 else 0.0
                
                pad_trend[year] = {
                    'infested_ha': infested_std_ha,
                    'infested_pct': infested_std_pct
                }
                print(f"    [+] Year {year}: {infested_std_ha:.2f} ha ({infested_std_pct:.1f}%)")
            else:
                pad_trend[year] = {'infested_ha': np.nan, 'infested_pct': np.nan}
                
        historical_trends[pad['name']] = pad_trend
        
        summary_row = {
            'Paddock_Name': pad['name'],
            'Nominal_Area_Ha': round(pad['nominal_area_ha'], 2),
            'Buffered_Area_Ha': round(total_buffered_ha, 2),
            'Targeted_Infested_Ha': round(infested_ha, 2),
            'Targeted_Infested_Pct': round(infestation_rate, 1),
            'Tier1_Eradication_Ha': round(t1_ha, 2) if not np.isnan(t1_ha) else np.nan,
            'Tier1_Eradication_Pct': round(t1_pct, 1) if not np.isnan(t1_pct) else np.nan,
            'Tier2_Containment_Ha': round(t2_ha, 2) if not np.isnan(t2_ha) else np.nan,
            'Tier2_Containment_Pct': round(t2_pct, 1) if not np.isnan(t2_pct) else np.nan,
            'Tier3_Refugia_Ha': round(t3_ha, 2) if not np.isnan(t3_ha) else np.nan,
            'Tier3_Refugia_Pct': round(t3_pct, 1) if not np.isnan(t3_pct) else np.nan,
        }
        # Add historical pct values directly to summary row
        for year in years:
            summary_row[f'Infested_Pct_{year}'] = pad_trend[year]['infested_pct']
            
        summary_data.append(summary_row)
        
    df_summary = pd.DataFrame(summary_data)
    summary_csv = os.path.join(out_dir, "paddock_infestation_summary.csv")
    df_summary.to_csv(summary_csv, index=False)
    print(f"\n[+] Saved paddock summary CSV to: {summary_csv}")
    
    # 3. Generate charts
    # Chart 1: Current targeted comparison
    plt.figure(figsize=(10, 6))
    colors = ['#ef4444' if pct > 40 else '#f97316' if pct > 15 else '#22c55e' for pct in df_summary['Targeted_Infested_Pct']]
    bars = plt.bar(df_summary['Paddock_Name'], df_summary['Targeted_Infested_Pct'], color=colors, edgecolor='#444', linewidth=1.2)
    plt.title('Targeted Buffel Grass Infestation Rate by Paddock', fontsize=14, fontweight='bold', pad=15)
    plt.ylabel('Infestation Rate (%)', fontsize=12)
    plt.xlabel('Paddock Name', fontsize=12)
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2.0, height + 1.5, f"{height:.1f}%", ha='center', fontsize=10, fontweight='bold')
        
    plt.ylim(0, max(df_summary['Targeted_Infested_Pct']) + 10)
    plt.xticks(rotation=15)
    plt.tight_layout()
    chart1_path = os.path.join(out_dir, "paddock_comparison_chart.png")
    plt.savefig(chart1_path, dpi=150)
    plt.close()
    
    # Chart 2: Historical trend per paddock
    plt.figure(figsize=(12, 7))
    line_styles = ['-', '--', '-.', ':', '-']
    markers = ['o', 's', '^', 'd', 'v']
    color_palette = ['#ef4444', '#3b82f6', '#10b981', '#f59e0b', '#8b5cf6']
    
    for idx, pad in enumerate(paddocks):
        name = pad['name']
        trend_y = [historical_trends[name][y]['infested_pct'] for y in years]
        plt.plot(years, trend_y, label=name, marker=markers[idx % len(markers)], 
                 linestyle=line_styles[idx % len(line_styles)], color=color_palette[idx % len(color_palette)],
                 linewidth=2.5, markersize=7)
                 
    plt.title('Historical Decade Infestation Trend (Paddock-by-Paddock)', fontsize=14, fontweight='bold', pad=15)
    plt.xlabel('Year', fontsize=12)
    plt.ylabel('Infestation Rate (%)', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xticks(years)
    plt.legend(fontsize=10, loc='upper left')
    plt.tight_layout()
    chart2_path = os.path.join(out_dir, "paddock_historical_trend.png")
    plt.savefig(chart2_path, dpi=150)
    plt.close()
    print(f"[+] Saved comparison charts to paddock_breakdown/ directory.")
    
    # 4. Generate Leaflet HTML Map with boundaries and popups
    leaflet_data = []
    for pad in paddocks:
        # Extract coordinates list
        poly_coords = []
        x_coords, y_coords = pad['poly_wgs84'].exterior.coords.xy
        for x, y in zip(x_coords, y_coords):
            poly_coords.append([y, x]) # [lat, lon]
            
        pad_row = df_summary[df_summary['Paddock_Name'] == pad['name']].iloc[0]
        
        # Color coding: Green (<15%), Orange (15-40%), Red (>40%)
        color = '#ef4444' if pad_row['Targeted_Infested_Pct'] > 40 else '#f97316' if pad_row['Targeted_Infested_Pct'] > 15 else '#22c55e'
        
        leaflet_data.append({
            'name': pad['name'],
            'coords': poly_coords,
            'color': color,
            'nominal_area': float(pad_row['Nominal_Area_Ha']),
            'buffered_area': float(pad_row['Buffered_Area_Ha']),
            'infested_area': float(pad_row['Targeted_Infested_Ha']),
            'infested_pct': float(pad_row['Targeted_Infested_Pct']),
            't1_ha': None if np.isnan(pad_row['Tier1_Eradication_Ha']) else float(pad_row['Tier1_Eradication_Ha']),
            't1_pct': None if np.isnan(pad_row['Tier1_Eradication_Pct']) else float(pad_row['Tier1_Eradication_Pct']),
            't2_ha': None if np.isnan(pad_row['Tier2_Containment_Ha']) else float(pad_row['Tier2_Containment_Ha']),
            't2_pct': None if np.isnan(pad_row['Tier2_Containment_Pct']) else float(pad_row['Tier2_Containment_Pct']),
            't3_ha': None if np.isnan(pad_row['Tier3_Refugia_Ha']) else float(pad_row['Tier3_Refugia_Ha']),
            't3_pct': None if np.isnan(pad_row['Tier3_Refugia_Pct']) else float(pad_row['Tier3_Refugia_Pct']),
            'trend': {int(y): float(historical_trends[pad['name']][y]['infested_pct']) for y in years if not np.isnan(historical_trends[pad['name']][y]['infested_pct'])}
        })
        
    # Get centers
    all_lats = []
    all_lons = []
    for pad in paddocks:
        all_lats.append(pad['centroid'][0])
        all_lons.append(pad['centroid'][1])
    center_lat = sum(all_lats) / len(all_lats)
    center_lon = sum(all_lons) / len(all_lons)
    
    # Load raster bounds for overlay
    da_gps = da_ref.rio.reproject("EPSG:4326")
    bounds = da_gps.rio.bounds()
    leaflet_bounds = [[bounds[1], bounds[0]], [bounds[3], bounds[2]]]
    
    # Render direct Leaflet HTML Map
    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <title>Keewong Property - Paddock breakdown Map</title>
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
            padding: 12px;
            background: rgba(30, 30, 36, 0.95);
            box-shadow: 0 4px 15px rgba(0,0,0,0.5);
            border-radius: 8px;
            border: 1px solid #444;
            color: #f3f3f3;
            max-width: 280px;
        }}
        .info-panel h2 {{
            margin: 0 0 8px 0;
            font-size: 14px;
            color: #38bdf8;
            border-bottom: 1px solid #444;
            padding-bottom: 4px;
        }}
        .paddock-popup {{
            font-size: 12px;
            color: #333;
            min-width: 220px;
        }}
        .paddock-popup h3 {{
            margin: 0 0 6px 0;
            font-size: 13px;
            color: #0369a1;
            border-bottom: 1px solid #ddd;
            padding-bottom: 3px;
        }}
        .paddock-popup table {{
            width: 100%;
            border-collapse: collapse;
        }}
        .paddock-popup td {{
            padding: 3px 0;
        }}
        .paddock-popup .val {{
            font-weight: bold;
            text-align: right;
        }}
        .paddock-popup .highlight {{
            color: #b91c1c;
            font-weight: bold;
        }}
        .legend {{
            padding: 8px;
            background: rgba(30, 30, 36, 0.95);
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

    var paddockData = {json.dumps(leaflet_data)};
    
    // Add overlays
    var maxSpikeOverlay = L.imageOverlay('../targeted/buffel_max_growth_spike.png', {leaflet_bounds}, {{
        opacity: 0.55,
        interactive: false
    }}).addTo(map);

    // Render paddocks
    paddockData.forEach(function(pad) {{
        var poly = L.polygon(pad.coords, {{
            color: pad.color,
            fillColor: pad.color,
            fillOpacity: 0.15,
            weight: 3
        }}).addTo(map);
        
        // Generate trend list html
        var trendHtml = "";
        Object.keys(pad.trend).sort().forEach(function(yr) {{
            trendHtml += yr + ": " + pad.trend[yr].toFixed(1) + "% | ";
        }});
        if (trendHtml.length > 3) {{
            trendHtml = trendHtml.substring(0, trendHtml.length - 3);
        }}

        var tierRows = "";
        if (pad.t1_ha !== null) {{
            tierRows = "<tr style='border-top:1px solid #eee;'><td>Tier 1 (Eradication):</td><td class='val' style='color:#16a34a;'>" + pad.t1_ha.toFixed(1) + " ha (" + pad.t1_pct.toFixed(1) + "%)</td></tr>" +
                "<tr><td>Tier 2 (Containment):</td><td class='val' style='color:#dc2626;'>" + pad.t2_ha.toFixed(1) + " ha (" + pad.t2_pct.toFixed(1) + "%)</td></tr>" +
                "<tr><td>Tier 3 (Refugia):</td><td class='val' style='color:#0891b2;'>" + pad.t3_ha.toFixed(1) + " ha (" + pad.t3_pct.toFixed(1) + "%)</td></tr>";
        }}

        var popupContent = "<div class='paddock-popup'>" +
            "<h3>" + pad.name + "</h3>" +
            "<table>" +
            "<tr><td>Nominal Area:</td><td class='val'>" + pad.nominal_area.toFixed(1) + " ha</td></tr>" +
            "<tr><td>Buffered Area:</td><td class='val'>" + pad.buffered_area.toFixed(1) + " ha</td></tr>" +
            "<tr><td>Core Infested:</td><td class='val highlight'>" + pad.infested_area.toFixed(1) + " ha (" + pad.infested_pct.toFixed(1) + "%)</td></tr>" +
            tierRows +
            "</table>" +
            "<div style='margin-top:8px; border-top:1px solid #ccc; padding-top:4px; font-size:10px; color:#555;'>" +
            "<strong>Historical Trend:</strong><br>" + trendHtml +
            "</div>" +
            "</div>";
            
        poly.bindPopup(popupContent);
    }});

    // Add legend
    var legend = L.control({{ position: 'bottomright' }});
    legend.onAdd = function (map) {{
        var div = L.DomUtil.create('div', 'legend');
        div.innerHTML = '<h3 style="margin:0 0 5px 0; font-size:12px; color:#38bdf8;">Infestation Level</h3>' +
            '<i style="background: #ef4444"></i> High (> 40%)<br>' +
            '<i style="background: #f97316"></i> Moderate (15 - 40%)<br>' +
            '<i style="background: #22c55e"></i> Low (< 15%)<br>';
        return div;
    }};
    legend.addTo(map);

    // Add info control
    var info = L.control({{ position: 'bottomleft' }});
    info.onAdd = function(map) {{
        var div = L.DomUtil.create('div', 'info-panel');
        div.innerHTML = '<h2>Paddock-by-Paddock Breakdown</h2>' +
            '<p style="font-size:11px; margin:0; color:#aaa;">Click on any paddock boundary line to view detailed infestation statistics and decade trend lines.</p>';
        return div;
    }};
    info.addTo(map);

</script>
</body>
</html>
"""
    
    html_path = os.path.join(out_dir, "paddock_breakdown_map.html")
    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    print(f"[+] Saved paddock Leaflet HTML map to: {html_path}")
    
    print("\n" + "=" * 60)
    print("[+] PADDOCK-BY-PADDOCK BREAKDOWN COMPLETE!")
    print(f"    Summary CSV: {summary_csv}")
    print(f"    Interactive Map: {html_path}")
    print("=" * 60)

if __name__ == "__main__":
    main()
