#!/bin/bash

set -x

THREADS=(16 32)

GRID_FILES=("icon_grid_0051_R02B05_G.nc")
GRID_NAMES=("R02B05")

# GRID_FILES=("icon_grid_0052_R02B06_G.nc")
# GRID_NAMES=("R02B06")

GRID_FILES=("icon_grid_0050_R02B03_G.nc" "icon_grid_0049_R02B04_G.nc")
GRID_NAMES=("R02B03" "R02B04")

# GRID_FILES=("icon_grid_0051_R02B05_G.nc" "icon_grid_0052_R02B06_G.nc")
# GRID_FILES=("icon_grid_0053_R02B07_G.nc")
# GRID_NAMES=("R02B05" "R02B06")
# GRID_NAMES=("R02B07")

for F2D_THREADS in "${THREADS[@]}"; do
  export F2D_THREADS
  for i in "${!GRID_FILES[@]}"; do
    F2D_GRID="/scratch/pmazumde/icon-grids/${GRID_FILES[$i]}"
    export F2D_GRID
    GRID_NAME="${GRID_NAMES[$i]}"

    echo "Running with F2D_THREADS=$F2D_THREADS, F2D_GRID=$F2D_GRID"
    bash exp.exclaim_ape_R02B04.run #>~/"${GRID_NAME}.t=${F2D_THREADS}.dace.log" 2>&1
  done
done

