---
title: "Readme"
format:
  html:
    embed-resources: true
---

## DisPerSE Workflow Cheatsheet

Concise reference for pulling Quijote data, running the DisPerSE pipelines, and exporting visualization-ready artifacts.

### Documentation
- Core pipelines: `docs/ANALYZE_SNAPSHOT_USER_GUIDE.qmd` (3D) and `docs/ANALYZE_SNAPSHOT_2D_USER_GUIDE.qmd` (2D projections)
- Slab + PNGs: `docs/BATCH_CLIP_USER_GUIDE.qmd`
- Batch tiling + per-slab: `docs/BATCH_CROP_AND_CLIP_USER_GUIDE.qmd`
- Topology overlap stats: `docs/ND_TOPO_STATS_USER_GUIDE.qmd`
- End-to-end runbook: `docs/TOPOLOGY_ANALYSIS_USER_GUIDE.qmd`
- Methods summary: `docs/PROJECT_METHODS.qmd`
- Utilities/other tools: `docs/OTHER_TOOLS_USER_GUIDE.qmd`
- Results examples: `docs/results.qmd` / `docs/results.html`

### Project Layout
- All scripts live in `scripts/`:
  - `scripts/analyze_snapshot.py` / `scripts/analyze_snapshot_2d.py`: main 3D/2D pipelines (decimate/crop, rebase for periodic Delaunay, MSE, convert to VTK/VTU/VTP with global shift-back).
  - `scripts/compute_density_field.py`: build a CIC density grid (HDF5).
  - `scripts/export_grid_vti.py`: wrap a 3D HDF5 grid as VTI for ParaView.
  - `scripts/export_snapshot_vtu.py`: stream snapshot particles to VTU (or VTI from an existing grid with `--density-field`).
  - `scripts/batch_clip.py`: clip slabs, flatten to 2D, average density, and render PNGs.
  - `scripts/ndtopo_stats.py`: compute overlap/unique/unassigned stats (uses Delaunay scalars, matches by `true_index`/`int(cell)`).
  - `scripts/batch_crop_and_clip.py`: tile a snapshot into many crops, run analyze → stats → slab clips for each.
- `data/`: place snapshots or density grids here. `outputs/`: per-run artifacts.

### Key Links
- DisPerSE docs: https://github.com/thierry-sousbie/DisPerSE | https://www.iap.fr/useriap/sousbie/web/html/indexd41d.html
- Quijote data: https://quijote-simulations.readthedocs.io/en/latest/access.html | https://quijote-simulations.readthedocs.io/en/latest/bsq.html
- Globus CLI: https://docs.globus.org/cli/

### Globus Transfers
Snapshot (BSQ/0):
```bash
globus endpoint local-id
export DST=d0b6b6f1-bdf0-11f0-9431-0e092d85c59b
globus endpoint search quijote
export SRC=e0eae0aa-5bca-11ea-9683-0e56c063f437
globus ls "$SRC:/Snapshots/fiducial/0/snapdir_000/"
globus transfer "$SRC:/Snapshots/fiducial/0/snapdir_000/" "$DST:~/Downloads/snapdir_000/" --label "Quijote Snapshot Fiducial 0 z0"
```
Snapshot (BSQ/0):
```bash
globus endpoint local-id
export DST=d0b6b6f1-bdf0-11f0-9431-0e092d85c59b
globus endpoint search quijote
export SRC=f4863854-3819-11eb-b171-0ee0d5d9299f
globus ls "$SRC:/Snapshots/BSQ/0/"
globus transfer "$SRC:/Snapshots/BSQ/0/snap_010.hdf5" "$DST:~/Downloads/snap_010.hdf5" --label "Quijote Snapshot BSQ 0 z0"
```
Density field (3D_cubes/BSQ/0):
```bash
globus endpoint local-id
export DST=d0b6b6f1-bdf0-11f0-9431-0e092d85c59b
globus endpoint search quijote
export SRC=e0eae0aa-5bca-11ea-9683-0e56c063f437
globus ls "$SRC:/3D_cubes/BSQ/0/"
globus transfer "$SRC:/3D_cubes/BSQ/0/df_m_CIC_z=0.00.hdf5" "$DST:~/Downloads/df_m_CIC_z=0.00.hdf5" --label "Quijote Density Field BSQ 0 z0"
```

### Snapshot → walls/filaments (DisPerSE)
Environment: `conda activate disperse` (with `hdf5plugin`, `vtkmodules`).

Typical run (cropped, periodic, rebased for tessellation, shifted VTK outputs):
```bash
python scripts/analyze_snapshot.py \
  --input data/snap_010.hdf5 \
  --output-dir outputs/snap_010_subbox \
  --output-prefix snap_010 \
  --crop-box 0 0 0 500000 500000 100000 \
  --stride 1 \ 
  --delaunay-btype periodic \
  --export-delaunay \
  --mse-nsig 3.0 \
  --dump-manifolds JE1a \
  --dump-arcs U \
  --netconv-smooth 20 \
  --skelconv-smooth 20
```
Notes:
- Coordinates are rebased to the crop origin for `-periodic` tessellation, then VTKs are shifted back to global.
- Outputs include shifted VTKs (`*_delaunay.vtu`, `*_manifolds_<TAG>.vtu`, `*_filaments_<ARC>.vtp`) plus per-file summary stats CSV.
- You can replace `--stride 1` with `--target-count 2000000` to auto-pick a stride near 2M particles.

### Density Workflows
From snapshot to density grid:
```bash
python scripts/compute_density_field.py \
  --input data/snap_010.hdf5 \
  --parttype PartType1 \
  --output outputs/snap_010_density.hdf5 \
  --grid-size 256 \
  --chunk-size 2000000 \
  --store-contrast
```
Visualize grid in ParaView (VTI):
```bash
python scripts/export_snapshot_vtu.py \
  --density-field \
  --input outputs/snap_010_density.hdf5 \
  --grid-dataset DensityField/density \
  --output outputs/snap_010_density.vti
```
Convert an existing grid directly:
```bash
python scripts/export_grid_vti.py \
  --input data/df_m_CIC_z_1_00.hdf5 \
  --dataset df \
  --output outputs/df_m_CIC_z_1_00.vti \
  --spacing 1 1 1 \
  --origin 0 0 0 \
  --field-name density
```

### Slab clip + PNGs
Clip a slab, flatten, average density, and render overlays:
```bash
PV_NO_MPI=1 /Applications/ParaView-6.0.1.app/Contents/bin/pvpython scripts/batch_clip.py \
  --input-dir outputs/snap_010_subbox \
  --walls snap_010_manifolds_JE1a.vtu \
  --filaments snap_010_filaments_U.vtp \
  --delaunay snap_010_delaunay.vtu \
  --output-dir outputs/snap_010_subbox/auto2d \
  --output-prefix snap_010 \
  --slab-axis z --slab-origin 0 --slab-thickness 10 \
  --resample-dims 512 512 64 \
  --scalar-name log_field_value \
  --save-pngs
```
Outputs: per-slab VTKs (3D + flattened + averaged density), PNGs, and summary stats.

### Topology overlap stats
Use Delaunay scalars to quantify shared/unique/unassigned topology membership:
```bash
python scripts/ndtopo_stats.py \
  --delaunay-vtk outputs/snap_010_subbox/snap_010_delaunay.vtu \
  --walls-vtk    outputs/snap_010_subbox/snap_010_manifolds_JE1a.vtu \
  --filaments-vtk outputs/snap_010_subbox/snap_010_filaments_U.vtp \
  --walls-id-field true_index \
  --filaments-id-field cell \
  --scalar-fields mass field_value log_field_value \
  --output-csv outputs/snap_010_subbox/snap_010_topology_stats.csv
```
`--write-vtk` will auto-convert NDnet/NDskl (unsmoothed) if VTKs are not present.

### Batch crops + slabs
Tile the box, run analyze → stats → slab clips for each crop:
```bash
python scripts/batch_crop_and_clip.py \
  --snapshot data/snap_010.hdf5 \
  --output-root outputs/quijote_batches \
  --crop-size 500000 500000 100000 \
  --x-range 0 1000000 --y-range 0 1000000 --z-range 0 200000 \
  --mse-nsig 5.0 \
  --dump-manifolds JE1a --dump-arcs U \
  --netconv-smooth 20 --skelconv-smooth 20 \
  --slab-step 10 --slab-thickness 10 \
  --resample-dims 500 500 100 \
  --scalar-name log_field_value
```
Outputs are organized per crop (`*_delaunay.vtu`, walls/filaments VTKs, per-crop topology stats) and per slab subfolder (flattened VTKs, PNGs, stats).

### Other Utilities
- CosmoFlow snapshot prep:
```bash
python scripts/prepare_cosmoflow_snapshot.py \
  --input data/univ_ics_2019-03_a10000668.hdf5 \
  --output data/cosmoflow_z0_snapshot.hdf5 \
  --redshift-index 0 \
  --target-count 1000000 \
  --threshold 5
```
- ParaView automation example:
```bash
/Applications/ParaView-6.0.1.app/Contents/bin/pvpython visualize_walls_paraview.py \
  --input outputs/snap_010/snap_010_walls.vtu \
  --field field_value \
  --threshold 0.001 0.02 \
  --scalar-bar-size 0.035 0.20 \
  --scalar-bar-position 0.85 0.1 \
  --scalar-bar-font-sizes 8 7 \
  --screenshot outputs/snap_010/walls_preview.png
```


## some AI queries (detailed guide)
> What is the benefit of DisPerSE? Could I just analyze the HDF5 in ParaView?

### Why DisPerSE vs. ParaView
- Reconstructs the Delaunay tessellation and Morse–Smale complex; persistence keeps only significant walls/filaments and honors periodic boxes.
- ParaView visualizes particles/grids but does no topology extraction; you’d only see raw density/points.
- `netconv`/`skelconv` convert DisPerSE outputs into VTK meshes you can inspect anywhere.

### Knobs that shape wall/filament geometry
- Particle sampling: `--stride` / `--target-count` (denser preserves fine structure; more thinning smooths).
- Boundaries: `--periodic`, `--delaunay-btype` (wraparound, block options).
+- Persistence: `--mse-nsig` or `--persistence-cut` (higher → prune weaker features; multiple values allowed, e.g., `--nsig 3.5 4.0 5.0`).
- Manifolds: `--dump-manifolds` (JD1d walls, J0a void minima, etc.).
- Threshold style: absolute (`--persistence-cut`) vs. sigma (`--mse-nsig`).
- Units: `--input-unit` / `--output-unit` to keep coordinate scales consistent.
- Post-visualization: `visualize_walls_paraview.py --threshold` can hide/show ranges.

If you’re experimenting, tweak in order: (1) input density/stride, (2) mse thresholds (-nsig vs -cut), (3) -dumpManifolds choice.

---

### Build/setup walkthrough (macOS)

**What we’re building**
- A reproducible setup that:
  - builds DisPerSE from source (CGAL, GMP/MPFR, GSL, CFITSIO),
  - downloads a Quijote snapshot (e.g., `snap_010.hdf5`),
  - optionally grids particles (CIC),
  - runs DisPerSE to extract filaments/walls,
  - converts results to VTK for ParaView.

**One-time prerequisites**
- Command Line Tools (if needed): `xcode-select --install`
- Homebrew libs (Apple Silicon paths assumed):
```bash
brew install cgal gmp mpfr gsl cfitsio boost@1.85
```
- Conda env (example):
```bash
conda create -n disperse -c conda-forge python=3.11 \
  numpy scipy h5py astropy tqdm numba psutil
conda activate disperse
# avoid HDF5 plugin crash:
unset HDF5_PLUGIN_PATH
export HDF5_USE_FILE_LOCKING=FALSE
```

**Build DisPerSE (once)**
```bash
cd /Users/fules/src/DisPerSE

# 1) Shim to ensure CGAL/GMP/MPFR are found
cat > cgal_use_dummy.cmake <<'CMAKE'
if (DEFINED CGAL_DIR AND EXISTS "${CGAL_DIR}/UseCGAL.cmake")
  include("${CGAL_DIR}/UseCGAL.cmake")
endif()
if (NOT GMP_FOUND)  find_package(GMP)  endif()
if (NOT MPFR_FOUND) find_package(MPFR) endif()
if (GMP_FOUND AND NOT TARGET GMP::gmp)
  add_library(GMP::gmp UNKNOWN IMPORTED)
  set_target_properties(GMP::gmp PROPERTIES
    IMPORTED_LOCATION "${GMP_LIBRARIES}"
    INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDE_DIR};${GMP_INCLUDE_DIRS}")
endif()
if (MPFR_FOUND AND NOT TARGET MPFR::mpfr)
  add_library(MPFR::mpfr UNKNOWN IMPORTED)
  set_target_properties(MPFR::mpfr PROPERTIES
    IMPORTED_LOCATION "${MPFR_LIBRARIES}"
    INTERFACE_INCLUDE_DIRECTORIES "${MPFR_INCLUDE_DIR};${MPFR_INCLUDE_DIRS}")
endif()
CMAKE

# 2) Configure & build
rm -rf build && mkdir build && cd build
cmake .. -Wno-dev \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
  -DGSL_DIR="$(brew --prefix gsl)" \
  -DCFITSIO_DIR="$(brew --prefix cfitsio)" \
  -DCGAL_DIR="$(brew --prefix cgal)/lib/cmake/CGAL" \
  -DBoost_NO_BOOST_CMAKE=ON \
  -DBoost_NO_SYSTEM_PATHS=ON \
  -DBOOST_ROOT="/opt/homebrew/opt/boost@1.85" \
  -DCMAKE_PREFIX_PATH="/opt/homebrew/opt/boost@1.85:/opt/homebrew" \
  -DCGAL_USE_FILE="$PWD/../cgal_use_dummy.cmake" \
  -DCMAKE_CXX_STANDARD=17 -DCMAKE_CXX_STANDARD_REQUIRED=ON -DCMAKE_CXX_EXTENSIONS=OFF \
  -DCMAKE_CXX_FLAGS="-D_LIBCPP_ENABLE_CXX17_REMOVED_UNARY_BINARY_FUNCTION -D_LIBCPP_ENABLE_CXX17_REMOVED_BINDERS"
make -j4  # builds: mse, skelconv, netconv, fieldconv, delaunay_2D/3D
```

If you see `Undefined symbols ___gmpn_*`, GMP/MPFR weren’t linked—this shim + Homebrew gmp/mpfr fixes it.

**PATH**
```bash
export PATH="/Users/fules/src/DisPerSE/build/src:$PATH"
```
Add to `~/.zshrc` if desired.

Docs: DisPerSE converters (skelconv) usage is on the DisPerSE site; the cosmology pipeline examples there mirror the `mse`/`skelconv`/`-manifolds` calls used here.

```bash
if (MPFR_FOUND AND NOT TARGET MPFR::mpfr)
  add_library(MPFR::mpfr UNKNOWN IMPORTED)
  set_target_properties(MPFR::mpfr PROPERTIES
    IMPORTED_LOCATION "${MPFR_LIBRARIES}"
    INTERFACE_INCLUDE_DIRECTORIES "${MPFR_INCLUDE_DIR};${MPFR_INCLUDE_DIRS}")
endif()
CMAKE

# 3.2 Configure & build
rm -rf build && mkdir build && cd build

cmake .. -Wno-dev \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
  -DGSL_DIR="$(brew --prefix gsl)" \
  -DCFITSIO_DIR="$(brew --prefix cfitsio)" \
  -DCGAL_DIR="$(brew --prefix cgal)/lib/cmake/CGAL" \
  -DBoost_NO_BOOST_CMAKE=ON \
  -DBoost_NO_SYSTEM_PATHS=ON \
  -DBOOST_ROOT="/opt/homebrew/opt/boost@1.85" \
  -DCMAKE_PREFIX_PATH="/opt/homebrew/opt/boost@1.85:/opt/homebrew" \
  -DCGAL_USE_FILE="$PWD/../cgal_use_dummy.cmake" \
  -DCMAKE_CXX_STANDARD=17 -DCMAKE_CXX_STANDARD_REQUIRED=ON -DCMAKE_CXX_EXTENSIONS=OFF \
  -DCMAKE_CXX_FLAGS="-D_LIBCPP_ENABLE_CXX17_REMOVED_UNARY_BINARY_FUNCTION -D_LIBCPP_ENABLE_CXX17_REMOVED_BINDERS"

make -j4  # builds: mse, skelconv, netconv, fieldconv, delaunay_2D/3D
```

If you ever see Undefined symbols ___gmpn_* again, it means CGAL needs GMP/MPFR. The shim above, plus Homebrew’s gmp/mpfr, fixes that.

Add the binaries to your shell PATH:

### Add to ~/.zshrc (or just export in the shell you use to run)
```bash
export PATH="/Users/fules/src/DisPerSE/build/src:$PATH"
```

Docs: DisPerSE converters (skelconv) usage is here (formats, flags).  ￼
A typical DisPerSE pipeline used in cosmological catalogs is illustrated here (you’ll recognize the mse/skelconv calls and the -manifolds flag).
