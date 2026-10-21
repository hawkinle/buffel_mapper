#!/usr/bin/env bash
# =============================================================================
# Standalone Linux Compiler for Buffel Mapper (PyInstaller)
# Produces: dist/buffel-mapper
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

VENV_DIR=".venv"
if [ ! -d "$VENV_DIR" ]; then
    echo "[*] Creating virtual environment in $VENV_DIR..."
    python3 -m venv "$VENV_DIR"
fi

echo "[*] Installing package and PyInstaller build dependencies..."
"$VENV_DIR/bin/pip" install --upgrade pip
"$VENV_DIR/bin/pip" install -e . pyinstaller

echo "[*] Compiling standalone Linux binary (dist/buffel-mapper)..."
"$VENV_DIR/bin/pyinstaller" \
    --noconfirm \
    --clean \
    --onefile \
    --name buffel-mapper \
    --paths src \
    --collect-all rasterio \
    --collect-all pyproj \
    --collect-all pyogrio \
    --collect-all geopandas \
    --collect-all rioxarray \
    --collect-all shapely \
    --collect-data xarray \
    --collect-data matplotlib \
    --hidden-import pystac_client \
    --hidden-import planetary_computer \
    --hidden-import PIL \
    src/buffel_mapper/cli.py

echo ""
echo "================================================================="
echo "[+] COMPILATION SUCCESSFUL!"
echo "    Standalone Linux Binary: $SCRIPT_DIR/dist/buffel-mapper"
echo "    Test with:               ./dist/buffel-mapper --help"
echo "================================================================="
