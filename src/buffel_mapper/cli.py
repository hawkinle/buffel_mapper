#!/usr/bin/env python3
"""
Unified CLI Entry Point for the Buffel Grass Remote Sensing Mapping Pipeline.
Designed for both standard Python execution and standalone PyInstaller compilation.
"""

import os
import sys
import shutil
import argparse

def _configure_bundled_proj_gdal():
    """Configure PROJ and GDAL data paths when running inside a compiled PyInstaller binary."""
    if getattr(sys, "frozen", False) and hasattr(sys, "_MEIPASS"):
        bundle_dir = sys._MEIPASS
        for candidate in [
            os.path.join(bundle_dir, "proj_data"),
            os.path.join(bundle_dir, "pyproj", "proj_dir", "share", "proj"),
            os.path.join(bundle_dir, "rasterio", "proj_data"),
        ]:
            if os.path.isdir(candidate):
                os.environ["PROJ_DATA"] = candidate
                os.environ["PROJ_LIB"] = candidate
                break
        for gdal_candidate in [
            os.path.join(bundle_dir, "gdal_data"),
            os.path.join(bundle_dir, "rasterio", "gdal_data"),
        ]:
            if os.path.isdir(gdal_candidate):
                os.environ["GDAL_DATA"] = gdal_candidate
                break

_configure_bundled_proj_gdal()

from buffel_mapper import (
    targeted_rain_analysis,
    generate_tier_map,
    historical_trend,
    spread_analyzer,
    paddock_breakdown,
    calibrate_model,
)


def copy_boundary_files(src_path, dst_dir):
    """Copy boundary KML or Shapefile (with all sidecar files) into the project directory."""
    ext = os.path.splitext(src_path)[1].lower()
    base_name = os.path.splitext(os.path.basename(src_path))[0]
    src_dir = os.path.dirname(os.path.abspath(src_path))

    if ext == ".kml":
        dst = os.path.join(dst_dir, "boundary.kml")
        shutil.copy(src_path, dst)
        return dst
    elif ext == ".shp":
        for sc in [".shp", ".shx", ".dbf", ".prj", ".cpg"]:
            sc_src = os.path.join(src_dir, base_name + sc)
            if os.path.exists(sc_src):
                shutil.copy(sc_src, os.path.join(dst_dir, "BSA" + sc))
        return os.path.join(dst_dir, "BSA.shp")
    else:
        dst = os.path.join(dst_dir, "boundary" + ext)
        shutil.copy(src_path, dst)
        return dst


def _invoke_module_main(module_main, argv_list):
    """Invoke a module's main() function in-process with the given CLI arguments."""
    old_argv = sys.argv[:]
    try:
        sys.argv = argv_list
        module_main()
    finally:
        sys.argv = old_argv


def cmd_run_all(args):
    """Run the complete end-to-end Buffel grass mapping workflow."""
    projects_root = os.path.abspath(args.projects_root)
    project_dir = os.path.join(projects_root, args.project)
    targeted_out = os.path.join(project_dir, "targeted")
    historical_out = os.path.join(project_dir, "historical")
    os.makedirs(targeted_out, exist_ok=True)
    os.makedirs(historical_out, exist_ok=True)

    print("=" * 65)
    print(f"[*] BUFFEL MAPPER WORKFLOW STARTING: {args.project}")
    print(f"[*] Boundary File:  {args.boundary}")
    print(f"[*] Rainfall File:  {args.rain}")
    print(f"[*] Project Output: {project_dir}")
    print("=" * 65)

    local_boundary = copy_boundary_files(args.boundary, project_dir)
    local_rain = os.path.join(project_dir, "rain_data.txt")
    shutil.copy(args.rain, local_rain)

    # Step 1: Targeted Rain-Pulse Analysis
    print("\n--- STEP 1: TARGETED WARM-SEASON RAIN-PULSE ANALYSIS ---")
    targeted_rain_analysis.run_targeted_analysis(
        kml_path=local_boundary,
        rain_path=local_rain,
        out_dir=targeted_out,
        threshold=args.threshold,
        max_cloud=args.max_cloud,
        index_name=args.index,
        buffer_pixels=args.buffer_pixels,
        warm_season_only=not args.include_winter,
    )

    # Step 2: 3-Tier Weed Management Delineation
    if os.path.exists(os.path.join(project_dir, "BSA.shp")):
        print("\n--- STEP 2: 3-TIER WEED MANAGEMENT DELINEATION ---")
        _invoke_module_main(
            generate_tier_map.main,
            [
                "generate_tier_map.py",
                "--project",
                project_dir,
                "--buffer-pixels",
                str(args.buffer_pixels),
            ],
        )

    # Step 3 & 4: Decade Historical Trend & Spatial Spread Vectors
    if not args.skip_historical:
        print("\n--- STEP 3: DECADE HISTORICAL VOLATILITY TREND ---")
        historical_trend.run_historical_trend(
            kml_path=local_boundary,
            start_year=args.start_year,
            end_year=args.end_year,
            max_cloud=args.max_cloud,
            out_dir=historical_out,
            threshold=0.12,
            buffer_pixels=args.buffer_pixels,
        )

        print("\n--- STEP 4: SPATIAL INVASION SPREAD VECTOR ANALYSIS ---")
        _invoke_module_main(
            spread_analyzer.main,
            [
                "spread_analyzer.py",
                "--kml",
                local_boundary,
                "--historical-dir",
                historical_out,
            ],
        )
        spread_png = os.path.join(historical_out, "infestation_spread_map.png")
        if os.path.exists(spread_png):
            shutil.copy(spread_png, os.path.join(project_dir, "infestation_spread_map.png"))

    # Step 5: Paddock-by-Paddock Breakdown
    if os.path.exists(os.path.join(project_dir, "BSA.shp")):
        print("\n--- STEP 5: PADDOCK-BY-PADDOCK BREAKDOWN ---")
        _invoke_module_main(
            paddock_breakdown.main,
            [
                "paddock_breakdown.py",
                "--project",
                project_dir,
                "--buffer-pixels",
                str(args.buffer_pixels),
            ],
        )

    print("\n" + "=" * 65)
    print(f"[+] WORKFLOW COMPLETE FOR PROJECT: {args.project}")
    print(f"    Project Directory: {project_dir}")
    print("=" * 65)


def main():
    parser = argparse.ArgumentParser(
        prog="buffel-mapper",
        description="Buffel Grass Remote Sensing Mapping & 3-Tier Management CLI",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # 1. `run` subcommand (End-to-End Pipeline)
    p_run = subparsers.add_parser("run", help="Run the complete end-to-end mapping workflow")
    p_run.add_argument("--boundary", "-b", required=True, help="Path to property boundary (.shp or .kml)")
    p_run.add_argument("--rain", "-r", required=True, help="Path to BoM rainfall log (.txt)")
    p_run.add_argument("--project", "-p", required=True, help="Project name (e.g. Keewong_2026)")
    p_run.add_argument("--projects-root", default="projects", help="Root directory for project outputs (default: ./projects)")
    p_run.add_argument("--index", default="msavi", choices=["msavi", "ndvi"], help="Vegetation index (default: msavi)")
    p_run.add_argument("--threshold", type=float, default=0.15, help="Rain-pulse growth spike threshold (default: 0.15)")
    p_run.add_argument("--max-cloud", type=float, default=10.0, help="Max cloud cover percentage (default: 10.0)")
    p_run.add_argument("--buffer-pixels", type=int, default=2, help="Boundary buffer in 10m pixels (default: 2)")
    p_run.add_argument("--include-winter", action="store_true", help="Include cool-season winter rain events (May-Sep)")
    p_run.add_argument("--skip-historical", action="store_true", help="Skip the 10-year annual historical volatility scan")
    p_run.add_argument("--start-year", type=int, default=2016, help="Historical trend start year (default: 2016)")
    p_run.add_argument("--end-year", type=int, default=2025, help="Historical trend end year (default: 2025)")

    # 2. `targeted` subcommand
    p_targ = subparsers.add_parser("targeted", help="Run Step 1: Targeted rain-pulse analysis")
    p_targ.add_argument("--kml", required=True, help="Path to boundary (.shp or .kml)")
    p_targ.add_argument("--rain-file", required=True, help="Path to BoM rainfall file")
    p_targ.add_argument("--out-dir", required=True, help="Output directory")
    p_targ.add_argument("--index", default="msavi", choices=["msavi", "ndvi"])
    p_targ.add_argument("--threshold", type=float, default=0.15)
    p_targ.add_argument("--max-cloud", type=float, default=10.0)
    p_targ.add_argument("--buffer-pixels", type=int, default=2)
    p_targ.add_argument("--include-winter", action="store_true")

    # 3. `tiers` subcommand
    p_tiers = subparsers.add_parser("tiers", help="Run Step 2: 3-Tier management layer generation")
    p_tiers.add_argument("--project", required=True, help="Project name or path to project directory")
    p_tiers.add_argument("--buffer-pixels", type=int, default=2)

    # 4. `paddocks` subcommand
    p_pad = subparsers.add_parser("paddocks", help="Run Step 5: Paddock-by-paddock breakdown")
    p_pad.add_argument("--project", required=True, help="Project name or path to project directory")
    p_pad.add_argument("--buffer-pixels", type=int, default=2)

    # 5. `calibrate` subcommand
    p_cal = subparsers.add_parser("calibrate", help="Run Step 6: Ground-truthing model calibration")
    p_cal.add_argument("--csv-path", required=True, help="Path to field validation points CSV")
    p_cal.add_argument("--project-dir", required=True, help="Path to project directory")
    p_cal.add_argument("--optimize", action="store_true", help="Optimize spectral thresholds against field points")

    args = parser.parse_args()

    if args.command == "run":
        cmd_run_all(args)
    elif args.command == "targeted":
        targeted_rain_analysis.run_targeted_analysis(
            args.kml,
            args.rain_file,
            args.out_dir,
            args.threshold,
            args.max_cloud,
            args.index,
            args.buffer_pixels,
            warm_season_only=not args.include_winter,
        )
    elif args.command == "tiers":
        _invoke_module_main(
            generate_tier_map.main,
            ["generate_tier_map.py", "--project", args.project, "--buffer-pixels", str(args.buffer_pixels)],
        )
    elif args.command == "paddocks":
        _invoke_module_main(
            paddock_breakdown.main,
            ["paddock_breakdown.py", "--project", args.project, "--buffer-pixels", str(args.buffer_pixels)],
        )
    elif args.command == "calibrate":
        cmd = ["calibrate_model.py", "-c", args.csv_path, "-p", args.project_dir]
        if args.optimize:
            cmd.append("--optimize")
        _invoke_module_main(calibrate_model.main, cmd)


if __name__ == "__main__":
    main()
