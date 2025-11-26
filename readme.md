
## Project layout
- `analyze_snapshot.py`: main pipeline that subsamples Gadget snapshots, runs the DisPerSE executables, and emits NDnet/VTU products.
- `visualize_walls_paraview.py`: ParaView helper script for loading the VTU walls and applying consistent coloring/thresholding.
- `scripts/`: misc utilities (e.g., `scripts/test.py`) or ad-hoc helpers.
- `data/`: raw Quijote snapshot inputs placed here before running the analysis.
- `outputs/`: per-snapshot directories containing NDfield/NDnet intermediates and final wall meshes.
- `__pycache__/`: Python bytecode cache created automatically by the interpreter.

## links:
https://github.com/thierry-sousbie/DisPerSE  
https://www.iap.fr/useriap/sousbie/web/html/indexd41d.html  

https://quijote-simulations.readthedocs.io/en/latest/access.html  
https://docs.globus.org/globus-connect-personal/install/mac/  
https://docs.globus.org/cli/  

https://quijote-simulations.readthedocs.io/en/latest/bsq.html  
https://quijote-simulations.readthedocs.io/en/latest/snapshots.html#snapshots  

## deal with globus transfer:
(disperse) fules@Fuless-MBP disperse % globus endpoint local-id  
d0b6b6f1-bdf0-11f0-9431-0e092d85c59b  
(disperse) fules@Fuless-MBP disperse % export DST=d0b6b6f1-bdf0-11f0-9431-0e092d85c59b  

(disperse) fules@Fuless-MBP disperse % globus endpoint search quijote  
(disperse) fules@Fuless-MBP disperse % export SRC=f4863854-3819-11eb-b171-0ee0d5d9299f  

(disperse) fules@Fuless-MBP disperse % globus ls "$SRC:/Snapshots/BSQ/0/"  
ICs/  
extra_files/  
Cosmo_params.dat  
done  
snap_006.hdf5  
snap_010.hdf5  

(disperse) fules@Fuless-MBP disperse % globus transfer \
  "$SRC:/Snapshots/BSQ/0/snap_010.hdf5" \
  "$DST:~/Downloads/snap_010.hdf5" \
  --label "Quijote Snapshot BSQ 0 z0"

## if not in disperse environment then set:
(base) fules@Fuless-MacBook-Pro disperse % conda activate disperse

## run the analysis
(disperse) fules@Fuless-MacBook-Pro disperse %      python analyze_snapshot.py \
         --input data/snap_010.hdf5 \
         --output-dir outputs/snap_010 \
         --periodic \
         --target-count 2000000 \
         --nsig 3.5 \
         --manifold-spec JD2d \
         --wall-format vtu

## visualize in paraview:
(disperse) fules@Fuless-MacBook-Pro disperse %     /Applications/ParaView-6.0.1.app/Contents/bin/pvpython visualize_walls_paraview.py \
    --input outputs/snap_010/snap_010_walls.vtu \
    --field field_value \
    --threshold 0.001 0.02 \
    --scalar-bar-size 0.035 0.20 \
    --scalar-bar-position 0.85 0.1 \
    --scalar-bar-font-sizes 8 7 \
    --screenshot outputs/snap_010/walls_preview.png

## some AI queries:
› what is the benefit of disperse in this workflow? could i analyze the hdf5 file directly in paraview?

• I’m explaining how DisPerSE differs from ParaView: it performs topological extraction like persistent homology and handles periodic boundaries, enabling detection of walls and filaments from particle data, whereas ParaView mainly provides visualization without topology analysis. I’ll walk through the pipeline steps and why DisPerSE’s approach is essential compared to manual thresholding or just loading snapshots in ParaView.

• DisPerSE is doing the heavy lifting here:

  - It takes the raw particle catalog and rebuilds the full Delaunay tessellation so it can estimate a continuous density field (via DTFE) without going onto a fixed grid.
  - From that tessellation it computes the Morse–Smale complex, applies persistence filtering (your --nsig 3.5), and extracts only the topologically significant descending manifolds of index‑1 saddles—i.e., coherent wall surfaces—respecting periodic boundaries.
  - netconv then converts those manifolds into an actual mesh file you can visualize anywhere.

  ParaView can open the HDF5 snapshot (or a converted VTK file) directly, but it will just show the particles or any scalar field you derive—there’s no in-built topology analysis to tell you which sheets/filaments are significant, stitch them into manifolds, or enforce persistence thresholds. You could, for example, grid the particles onto a volume and use isosurfaces/volume rendering, but you’d be looking at raw density, not the Morse–Smale walls that DisPerSE extracts. In short: ParaView handles visualization; DisPerSE supplies the scientifically meaningful structures to visualize.

  Key DisPerSE Inputs That Shape Wall Geometry

Particle sampling (--stride/--target-count) – lowering the stride gives DisPerSE a denser point set, which captures more small-scale walls; raising it smooths over fine structure. This is the first lever on the geometry because everything downstream works on the decimated catalog.

Smoothing scale (delaunay_3D options) – DisPerSE’s triangulation doesn’t have many tunables beyond periodic boundaries and block size, but enabling --periodic (which adds the -periodic flag) makes features wrap through the domain instead of breaking at edges, crucial for large cosmological boxes.

Persistence filtering (mse -nsig/-cut) – the main control described in the linked docs: higher --mse-nsig (or higher --persistence-cut values) prunes low-significance pairs, leaving only strong walls; lower thresholds keep more tenuous structures. You can supply one value per critical pair (e.g., --nsig 3.5 4.0 5.0) exactly as DisPerSE expects.

Manifold descriptor (--manifold-spec) – DisPerSE lets you pick what descending manifolds to dump (e.g., JD1d merges all index-1 saddles into wall sheets, BC_p2 might restrict to certain boundary conditions). Choosing a different descriptor changes which segments of the Morse-Smale complex end up in your mesh, per the “Manifolds” section on the website.

Absolute vs. relative filtering – using --persistence-cut (absolute density thresholds) rather than --nsig (sigma-based) ties wall retention directly to the scalar field’s units. That’s helpful if you know the physical density difference you care about.

Unit conversion (--input-unit/--output-unit) – while not a DisPerSE parameter per se, scaling coordinates affects the reported dimensions in the final NDnet/VTU; DisPerSE assumes coordinates come in the same units you specify, so ensuring you export in the desired comoving length gives walls the correct parsec/megaparsec extent.

Thresholding post-processing (visualize_walls_paraview.py --threshold) – this happens after DisPerSE but still affects perceived structure by hiding low/high values in the chosen scalar field, effectively mimicking an additional persistence cut during visualization.

If you want to experiment like the documentation suggests, focus on: (1) density of the input point cloud, (2) mse thresholds (-nsig versus -cut), and (3) the -dumpManifolds spec. Those are the knobs that most directly translate into thicker/thinner, sparser/denser wall networks when compared against the DisPerSE settings described at https://www.iap.fr/useriap/sousbie/web/html/indexd41d.html.



  
