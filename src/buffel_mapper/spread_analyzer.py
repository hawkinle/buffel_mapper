#!/usr/bin/env python3
import os
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image

def load_boundary_coordinates_leaflet(file_path):
    import geopandas as gpd
    import shapely.geometry
    
    ext = os.path.splitext(file_path)[1].lower()
    if ext == ".kml":
        tree = ET.parse(file_path)
        root = tree.getroot()
        namespaces = {
            'kml': 'http://www.opengis.net/kml/2.2',
            'gx': 'http://www.google.com/kml/ext/2.2'
        }
        coordinates_elements = root.findall('.//kml:coordinates', namespaces)
        if not coordinates_elements:
            coordinates_elements = root.findall('.//coordinates')
        coord_text = coordinates_elements[0].text.strip()
        coords = []
        for pt in coord_text.split():
            parts = pt.split(',')
            if len(parts) >= 2:
                lon, lat = float(parts[0]), float(parts[1])
                coords.append([lat, lon])
        return coords
    elif ext in [".shp", ".geojson", ".gpkg"]:
        gdf = gpd.read_file(file_path)
        if gdf.empty:
            raise ValueError("Vector file is empty.")
        gdf_wgs84 = gdf.to_crs("EPSG:4326")
        geom = gdf_wgs84.geometry.iloc[0]
        
        if isinstance(geom, shapely.geometry.MultiPolygon):
            coords = []
            for poly in geom.geoms:
                poly_coords = []
                x_coords, y_coords = poly.exterior.coords.xy
                for x, y in zip(x_coords, y_coords):
                    poly_coords.append([y, x])
                coords.append(poly_coords)
            return coords
        else:
            x_coords, y_coords = geom.exterior.coords.xy
            coords = []
            for x, y in zip(x_coords, y_coords):
                coords.append([y, x]) # [latitude, longitude]
            return coords
    else:
        raise ValueError(f"Unsupported boundary format: {ext}")

def haversine_distance(lat1, lon1, lat2, lon2):
    """
    Calculate the great circle distance between two points on the earth (in meters).
    """
    R = 6371000.0  # Earth radius in meters
    phi1 = np.radians(lat1)
    phi2 = np.radians(lat2)
    delta_phi = np.radians(lat2 - lat1)
    delta_lambda = np.radians(lon2 - lon1)
    
    a = np.sin(delta_phi/2.0)**2 + np.cos(phi1) * np.cos(phi2) * np.sin(delta_lambda/2.0)**2
    c = 2.0 * np.arctan2(np.sqrt(a), np.sqrt(1.0 - a))
    return R * c

def calculate_bearing(lat1, lon1, lat2, lon2):
    """
    Calculate the bearing (angle in degrees) between two points.
    """
    y = np.sin(np.radians(lon2 - lon1)) * np.cos(np.radians(lat2))
    x = np.cos(np.radians(lat1)) * np.sin(np.radians(lat2)) - \
        np.sin(np.radians(lat1)) * np.cos(np.radians(lat2)) * np.cos(np.radians(lon2 - lon1))
    bearing = np.degrees(np.arctan2(y, x))
    return (bearing + 360.0) % 360.0

def get_cardinal_direction(bearing):
    directions = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
                  "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"]
    index = int((bearing + 11.25) / 22.5) % 16
    return directions[index]

import argparse

def main():
    parser = argparse.ArgumentParser(description="Buffel Infestation Spread Centroid Vector Analyzer")
    parser.add_argument("--kml", type=str, default="/home/leith/Documents/test.kml", help="Path to KML boundary file")
    parser.add_argument("--historical-dir", type=str, default="/home/leith/AI/Buffel_mapping/output_historical", help="Path to historical trend directory")
    
    args = parser.parse_args()
    kml_path = args.kml
    historical_dir = args.historical_dir
    csv_path = os.path.join(historical_dir, "historical_infestation_trend.csv")
    
    if not os.path.exists(csv_path):
        print(f"[!] Historical trend data not found at: {csv_path}")
        return
        
    df = pd.read_csv(csv_path)
    results = df.to_dict('records')
    
    print("[*] Parsing property boundary...")
    leaflet_coords = load_boundary_coordinates_leaflet(kml_path)
    
    # Check if MultiPolygon (list of lists of lists)
    is_multipolygon = isinstance(leaflet_coords[0][0], list)
    
    if is_multipolygon:
        all_lats = []
        all_lons = []
        for poly_coords in leaflet_coords:
            all_lats.extend([c[0] for c in poly_coords])
            all_lons.extend([c[1] for c in poly_coords])
        min_lat, max_lat = min(all_lats), max(all_lats)
        min_lon, max_lon = min(all_lons), max(all_lons)
    else:
        lats = [c[0] for c in leaflet_coords]
        lons = [c[1] for c in leaflet_coords]
        min_lat, max_lat = min(lats), max(lats)
        min_lon, max_lon = min(lons), max(lons)
    
    # We will compute the spatial centroid of infestation for each year
    centroids = []
    
    print("[*] Processing annual maps to extract infestation centroids...")
    for r in results:
        year = int(r['Year'])
        img_name = r['Image_File']
        img_path = os.path.join(historical_dir, img_name)
        
        if not os.path.exists(img_path):
            print(f"      [-] Warning: Map image not found for year {year}: {img_path}")
            continue
            
        # Load PNG and extract the green channel as a proxy for the viridis colormap value
        img = Image.open(img_path)
        img_arr = np.array(img)
        
        # In viridis, yellow/green represents the higher NDVI volatility values (infested)
        # Green channel > 140 works well to isolate values >= 0.12 (infested)
        # We also need to ignore the background (which is black/transparent outside KML)
        green_channel = img_arr[..., 1]
        alpha_channel = img_arr[..., 3] if img_arr.shape[-1] == 4 else np.ones_like(green_channel) * 255
        
        infested_pixels = (green_channel > 130) & (alpha_channel > 50)
        
        # Calculate pixel coordinates of infested mask
        y_indices, x_indices = np.where(infested_pixels)
        
        if len(y_indices) == 0:
            print(f"      [-] Year {year}: No active infestation pixels detected (dormant year).")
            continue
            
        # Map pixel coordinates to geographic coordinates using simple linear interpolation
        h, w = green_channel.shape
        
        # Note: raster rows (y) go from top (max_lat) to bottom (min_lat)
        # raster columns (x) go from left (min_lon) to right (max_lon)
        pixel_lats = max_lat - (y_indices / h) * (max_lat - min_lat)
        pixel_lons = min_lon + (x_indices / w) * (max_lon - min_lon)
        
        centroid_lat = np.mean(pixel_lats)
        centroid_lon = np.mean(pixel_lons)
        
        centroids.append({
            'Year': year,
            'Lat': centroid_lat,
            'Lon': centroid_lon,
            'Area_Ha': r['Infested_Area_Ha']
        })
        print(f"      [+] Year {year} Centroid: ({centroid_lat:.5f}, {centroid_lon:.5f}) - Area: {r['Infested_Area_Ha']:.1f} ha")
        
    df_centroids = pd.DataFrame(centroids)
    
    # Calculate year-over-year movements
    movements = []
    print("\n[*] Calculating year-over-year movement vectors...")
    for i in range(len(df_centroids) - 1):
        row1 = df_centroids.iloc[i]
        row2 = df_centroids.iloc[i+1]
        
        y1, y2 = int(row1['Year']), int(row2['Year'])
        dist = haversine_distance(row1['Lat'], row1['Lon'], row2['Lat'], row2['Lon'])
        bearing = calculate_bearing(row1['Lat'], row1['Lon'], row2['Lat'], row2['Lon'])
        direction = get_cardinal_direction(bearing)
        
        # Average velocity in meters per year
        velocity = dist / (y2 - y1)
        
        movements.append({
            'From_Year': y1,
            'To_Year': y2,
            'Distance_m': dist,
            'Velocity_m_yr': velocity,
            'Bearing_deg': bearing,
            'Direction': direction
        })
        print(f"      [+] {y1} -> {y2}: Moved {dist:.1f} meters {direction} (Speed: {velocity:.1f} m/year)")
        
    df_movements = pd.DataFrame(movements)
    df_movements.to_csv(os.path.join(historical_dir, "infestation_movement_vectors.csv"), index=False)
    
    # Calculate average movement metrics (excluding years where infestation dropped to near-zero dormancy)
    valid_moves = df_movements[(df_movements['Distance_m'] > 20) & (df_movements['Velocity_m_yr'] < 2000)]
    if not valid_moves.empty:
        avg_vel = valid_moves['Velocity_m_yr'].mean()
        # Find net movement from first year with active infestation to last year
        active_years = df_centroids[df_centroids['Area_Ha'] > 100]
        if len(active_years) >= 2:
            first_row = active_years.iloc[0]
            last_row = active_years.iloc[-1]
            net_dist = haversine_distance(first_row['Lat'], first_row['Lon'], last_row['Lat'], last_row['Lon'])
            net_bearing = calculate_bearing(first_row['Lat'], first_row['Lon'], last_row['Lat'], last_row['Lon'])
            net_dir = get_cardinal_direction(net_bearing)
            print(f"\n[+] Net Active Invasion Spread: Moved {net_dist:.1f} meters {net_dir} from {int(first_row['Year'])} to {int(last_row['Year'])}")
            print(f"[+] Average Expansion Velocity: {avg_vel:.1f} meters per year during active seasons.")
            
    # Generate 2-Panel Plot (Map on Left, Bar Chart on Right)
    print("\n[*] Exporting side-by-side infestation spread dashboard plot...")
    fig, (ax_map, ax_chart) = plt.subplots(1, 2, figsize=(18, 8.5), facecolor='#0f172a')
    
    # ------------------ LEFT PANEL: MAP ------------------
    ax_map.set_facecolor('#0f172a')
    
    # Draw property boundary and fill paddocks
    paddocks = []
    if is_multipolygon:
        first = True
        for idx, poly_coords in enumerate(leaflet_coords):
            lats_p = [c[0] for c in poly_coords]
            lons_p = [c[1] for c in poly_coords]
            centroid_lat = sum(lats_p) / len(lats_p)
            centroid_lon = sum(lons_p) / len(lons_p)
            
            # Name geographically (consistent with paddock_breakdown.py)
            if centroid_lat > -30.570 and centroid_lon > 147.850:
                name = "Central-East"
            elif centroid_lat < -30.590 and centroid_lon > 147.865:
                name = "South-East"
            elif centroid_lat < -30.580 and centroid_lon > 147.850:
                name = "South-Central"
            elif centroid_lat < -30.570 and centroid_lon < 147.840:
                name = "Western"
            else:
                name = "North-Western"
                
            # Fill paddock with very subtle dark slate color
            ax_map.fill(lons_p, lats_p, color='#1e293b', alpha=0.8, zorder=1)
            # Outline paddock
            ax_map.plot(lons_p + [lons_p[0]], lats_p + [lats_p[0]], color='#475569', linewidth=1.5, zorder=2,
                        label='Paddock Boundary' if first else "")
            first = False
            
            # Label paddock name at its centroid
            ax_map.text(centroid_lon, centroid_lat, name, 
                        color='#64748b', fontsize=8, fontweight='bold', fontstyle='italic', 
                        ha='center', va='center', zorder=2, alpha=0.8)
    else:
        ax_map.fill(lons, lats, color='#1e293b', alpha=0.8, zorder=1)
        ax_map.plot(lons + [lons[0]], lats + [lats[0]], color='#475569', linewidth=2.0, label='Property Boundary', zorder=2)
        
    # Draw chronological trajectory dotted line connecting all centroids in sequence
    lons_path = df_centroids['Lon'].tolist()
    lats_path = df_centroids['Lat'].tolist()
    ax_map.plot(lons_path, lats_path, color='#94a3b8', linestyle='--', linewidth=1.8, zorder=3, alpha=0.7)
    
    # Plot centroids with constant small size to prevent overlapping size clutter
    scatter = ax_map.scatter(
        df_centroids['Lon'], df_centroids['Lat'],
        c=df_centroids['Year'], cmap='viridis',
        s=80, 
        edgecolors='#f8fafc', linewidths=1.2, alpha=0.9, zorder=5
    )
    
    # Annotate years directly next to each centroid point
    for idx, row in df_centroids.iterrows():
        y_offset = 0.0007 if idx % 2 == 0 else -0.0013
        x_offset = 0.0006
        ax_map.text(row['Lon'] + x_offset, row['Lat'] + y_offset, f"{int(row['Year'])}", 
                    color='#f8fafc', fontsize=8.5, fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='#0f172a', alpha=0.75, edgecolor='none'),
                    zorder=10)
                 
    # Draw Net Spread Vector arrow highlighting the overall expansion trend
    c_first = df_centroids.iloc[0]
    c_last = df_centroids.iloc[-1]
    ax_map.annotate(
        f"Net Spread ({int(c_first['Year'])}->{int(c_last['Year'])})",
        xy=(c_last['Lon'], c_last['Lat']),
        xytext=(c_first['Lon'], c_first['Lat']),
        arrowprops=dict(arrowstyle="-|>", color='#ef4444', lw=2.5, mutation_scale=15, shrinkA=5, shrinkB=5),
        color='#f87171', fontsize=8.5, fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.2', facecolor='#0f172a', alpha=0.85, edgecolor='none'),
        zorder=4
    )
    
    # Add a North Arrow in the upper-right corner
    arrow_x = max_lon - 0.005
    arrow_y = max_lat - 0.002
    ax_map.annotate('N', xy=(arrow_x, arrow_y), xytext=(arrow_x, arrow_y - 0.0025),
                    arrowprops=dict(facecolor='#f8fafc', width=2, headwidth=6, headlength=6),
                    ha='center', va='center', color='#f8fafc', fontsize=8.5, fontweight='bold', zorder=10)
                 
    # Calculate scale bar length dynamically for 1 km
    mean_lat = (min_lat + max_lat) / 2.0
    meters_per_degree_lon = 111320.0 * np.cos(np.radians(mean_lat))
    scale_bar_len_deg = 1000.0 / meters_per_degree_lon
    
    # Place scale bar in the bottom-left corner
    scale_x_start = min_lon - 0.004
    scale_x_end = scale_x_start + scale_bar_len_deg
    scale_y = min_lat - 0.004
    
    ax_map.plot([scale_x_start, scale_x_end], [scale_y, scale_y], color='#f8fafc', linewidth=2.5, zorder=10)
    ax_map.text((scale_x_start + scale_x_end)/2.0, scale_y + 0.0006, "1 km", 
                color='#f8fafc', fontsize=8, ha='center', fontweight='bold', zorder=10)
             
    # Set Map bounds with balanced padding
    ax_map.set_xlim(min_lon - 0.008, max_lon + 0.015)
    ax_map.set_ylim(min_lat - 0.008, max_lat + 0.006)
    
    ax_map.set_title('A. Spatial Centroid Spread Trajectory', fontsize=12, fontweight='bold', pad=12, color='#f8fafc')
    ax_map.set_xlabel('Longitude', fontsize=9, color='#64748b')
    ax_map.set_ylabel('Latitude', fontsize=9, color='#64748b')
    
    # Style axis lines (spines) and grid
    for spine in ax_map.spines.values():
        spine.set_color('#334155')
    ax_map.tick_params(colors='#64748b', labelsize=8)
    ax_map.grid(True, linestyle=':', alpha=0.15)
    
    # ------------------ RIGHT PANEL: BAR CHART ------------------
    ax_chart.set_facecolor('#0f172a')
    
    years_str = [str(int(y)) for y in df_centroids['Year']]
    areas_ha = df_centroids['Area_Ha'].tolist()
    
    # Generate matching colors from the Viridis colormap dynamically
    cmap = plt.get_cmap('viridis')
    norm = plt.Normalize(df_centroids['Year'].min(), df_centroids['Year'].max())
    bar_colors = [cmap(norm(y)) for y in df_centroids['Year']]
    
    bars = ax_chart.bar(years_str, areas_ha, color=bar_colors, edgecolor='#334155', width=0.6, zorder=3)
    
    # Add values on top of bars
    for bar in bars:
        height = bar.get_height()
        ax_chart.text(
            bar.get_x() + bar.get_width()/2.0, 
            height + (max(areas_ha)*0.015),
            f"{height:.0f} ha",
            ha='center', va='bottom', fontsize=8, color='#f8fafc', fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.15', facecolor='#1e293b', alpha=0.8, edgecolor='none')
        )
        
    ax_chart.set_title('B. Annual Infestation Magnitude', fontsize=12, fontweight='bold', pad=12, color='#f8fafc')
    ax_chart.set_xlabel('Year of Analysis', fontsize=9, color='#64748b')
    ax_chart.set_ylabel('Estimated Infested Area (Hectares)', fontsize=9, color='#64748b')
    
    # Style spines and ticks
    for spine in ax_chart.spines.values():
        spine.set_color('#334155')
    ax_chart.tick_params(colors='#64748b', labelsize=8)
    ax_chart.grid(True, axis='y', linestyle=':', alpha=0.15)
    ax_chart.set_ylim(0, max(areas_ha) * 1.1)
    
    # Set main figure title
    plt.suptitle('Keewong Property - Buffel Grass Infestation Dynamics Dashboard (2016-2025)', 
                 fontsize=14, fontweight='bold', color='#f8fafc', y=0.97)
    
    plot_path = os.path.join(historical_dir, "infestation_spread_map.png")
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    plt.savefig(plot_path, dpi=150, facecolor='#0f172a')
    plt.close()
    print(f"[+] Saved spread vector map: {plot_path}")

if __name__ == "__main__":
    main()
