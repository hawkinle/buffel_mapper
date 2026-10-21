# 03. Compilation & Deployment Guide (Linux)

**Navigation:** [Home](Home.md) | [01. Methodology](01-Ecology-and-Remote-Sensing-Methodology.md) | [02. Architecture & CLI](02-Pipeline-Architecture-and-CLI-Reference.md) | **03. Compilation** | [04. Field Validation](04-Field-Validation-and-Calibration-Protocol.md)

---

## 1. Overview of the Compilation Process

Geospatial Python applications rely on native C/C++ shared libraries (`libgdal`, `libproj`, `libgeos`) and coordinate reference system databases (`proj.db`). Standard Python bundlers often fail at runtime unless these binary assets and environment paths are explicitly packaged.

`buffel-mapper` includes a turnkey Linux compilation pipeline:
1. **`src/buffel_mapper/cli.py` (`_configure_bundled_proj_gdal`):** Automatically detects when the application is executing inside a frozen PyInstaller bundle (`sys._MEIPASS`) and points `PROJ_DATA`, `PROJ_LIB`, and `GDAL_DATA` to the internal bundled database at startup.
2. **`build_executable.sh`:** Collects all shared libraries and data files from `rasterio`, `pyproj`, `pyogrio`, `geopandas`, `rioxarray`, `shapely`, `xarray`, and `matplotlib` into a single standalone Linux executable (`dist/buffel-mapper`).

---

## 2. How to Compile the Standalone Linux Binary

From the root of the repository, run:

```bash
chmod +x build_executable.sh
./build_executable.sh
```

### What `build_executable.sh` Does:
1. Ensures `.venv` exists and installs `buffel-mapper` alongside `pyinstaller`.
2. Runs `pyinstaller --onefile` targeting `src/buffel_mapper/cli.py`.
3. Outputs the standalone executable to:
   ```text
   dist/buffel-mapper
   ```

---

## 3. Running & Verifying the Compiled Executable

Verify that the compiled binary loads all geospatial libraries cleanly:

```bash
./dist/buffel-mapper --help
```

Run the full workflow using the compiled binary from any directory:

```bash
./dist/buffel-mapper run \
  --boundary /path/to/BSA.shp \
  --rain /path/to/rain_data.txt \
  --project Keewong_Compiled_Run
```
