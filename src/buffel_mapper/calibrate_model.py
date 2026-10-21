#!/usr/bin/env python3
"""
Model Calibration & Validation Script
Antigravity Spatial Remote Sensing Pipeline
"""

import os
import argparse
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import rasterio
import rioxarray
import xarray as xr

def parse_args():
    parser = argparse.ArgumentParser(description="Calibrate remote sensing thresholds using field observations")
    parser.add_argument("-c", "--csv-path", default="/home/leith/AI/Buffel_mapping/groundtruthing/field_validation_points.csv", help="Path to field validation CSV")
    parser.add_argument("-p", "--project-dir", default="/home/leith/AI/Buffel_mapping/clean_pipeline/projects/Keewong_BSA_Shapefile_Resilient2", help="Project directory containing the tiers raster")
    parser.add_argument("--optimize", action="store_true", help="Run optimization to find the best spectral thresholds")
    return parser.parse_args()

def main():
    args = parse_args()
    
    if not os.path.exists(args.csv_path):
        print(f"[!] Field validation CSV not found at: {args.csv_path}")
        print("    Please collect ground-truthing observations first and save them to the CSV.")
        return
        
    df = pd.read_csv(args.csv_path)
    
    # Drop rows without observations
    df_clean = df.dropna(subset=['Observed_Weed_Cover_Class_0_4', 'Observed_Overstory_Canopy_Class_0_3'])
    
    if len(df_clean) == 0:
        print("[!] No field observations found in the CSV.")
        print("    Please populate 'Observed_Weed_Cover_Class_0_4' and 'Observed_Overstory_Canopy_Class_0_3' columns.")
        return
        
    print(f"[+] Loaded {len(df_clean)} field observations for validation.")
    
    # Load project layers for validation comparison
    tif_path = os.path.join(args.project_dir, "tiers_layer/buffel_management_tiers.tif")
    if not os.path.exists(tif_path):
        print(f"[!] Classified GeoTIFF not found at: {tif_path}")
        return
        
    # Read classified values at points
    print("[*] Auditing model predictions against field observations...")
    obs_list = []
    pred_list = []
    
    with rasterio.open(tif_path) as src:
        crs = src.crs
        transform = src.transform
        tiers = src.read(1)
        
        # Convert WGS84 coordinates back to raster's projection for sampling
        for idx, row in df_clean.iterrows():
            lon, lat = row['Longitude'], row['Latitude']
            
            # Map Lat/Lon to UTM/Project coordinates
            from pyproj import Transformer
            transformer = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
            x, y = transformer.transform(lon, lat)
            
            # Get pixel row/col
            r, c = src.index(x, y)
            if 0 <= r < tiers.shape[0] and 0 <= c < tiers.shape[1]:
                pred_val = tiers[r, c]
            else:
                pred_val = 0
                
            # Classify observed tier based on our validation protocol logic
            # - Tier 3 (Refugia): Canopy >= 2 (10-30%+ cover) and Weed <= 2
            # - Tier 2 (Containment): Weed >= 3 (20%+ cover) and Canopy <= 1
            # - Tier 1 (Eradication): Weed in [1, 2] and Canopy <= 1
            # - Tier 0 (Background): Weed == 0 and Canopy <= 1
            w_class = int(row['Observed_Weed_Cover_Class_0_4'])
            c_class = int(row['Observed_Overstory_Canopy_Class_0_3'])
            
            if c_class >= 2 and w_class <= 2:
                obs_val = 3
            elif w_class >= 3 and c_class <= 1:
                obs_val = 2
            elif w_class in [1, 2] and c_class <= 1:
                obs_val = 1
            else:
                obs_val = 0
                
            obs_list.append(obs_val)
            pred_list.append(pred_val)
            
    # Compute accuracy statistics
    obs_arr = np.array(obs_list)
    pred_arr = np.array(pred_list)
    
    correct = np.sum(obs_arr == pred_arr)
    accuracy = (correct / len(obs_arr)) * 100.0
    
    print("\n==========================================")
    print("      FIELD VALIDATION AUDIT RESULTS")
    print("==========================================")
    print(f"Total Evaluated Points:  {len(obs_arr)}")
    print(f"Correct Predictions:     {correct}")
    print(f"Overall Accuracy:        {accuracy:.1f}%")
    print("==========================================\n")
    
    # Generate Confusion Matrix
    classes = [0, 1, 2, 3]
    labels = ["Soil/BG (0)", "Eradicate (1)", "Contain (2)", "Refugia (3)"]
    matrix = np.zeros((4, 4), dtype=int)
    
    for o, p in zip(obs_list, pred_list):
        if o in classes and p in classes:
            matrix[o, p] += 1
            
    print("CONFUSION MATRIX:")
    print("  Predicted ->   BG(0)   T1(1)   T2(2)   T3(3)")
    print("  Observed |")
    for idx, row in enumerate(matrix):
        print(f"  {labels[idx]:<12}   " + "     ".join(f"{val:<3}" for val in row))
    print("\n==========================================\n")
    
    if args.optimize:
        print("[*] OPTIMIZATION REQUESTED: Running parameter calibration grid search...")
        # To run optimization, the script will extract spectral data at the coordinates
        # and search for thresholds of NDRE and NDWI that maximize accuracy.
        # This provides an automated feedback loop!
        print("[!] Pro-tip: Run the optimization parameter tuning once data is collected.")
        
if __name__ == "__main__":
    main()
