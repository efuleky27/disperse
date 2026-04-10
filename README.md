# Disentangling the Cosmic Web

This repository contains my work on the project that I will present at the 2026 Regeneron International Science and Engineering Fair in Phoenix, Arizona.
I created an automated pipeline for extracting, analyzing, and visualizing the cosmic web — walls, filaments, and clusters — from N-body cosmological simulations.

## Documentation

Full documentation, workflow guides, and AI session logs are available at:

**<https://efuleky27.github.io/disperse/>**

## Overview

The pipeline takes Gadget/Quijote HDF5 snapshots, runs Morse–Smale topology extraction via DisPerSE, produces per-crop statistics and 2D/3D visualizations, and aggregates results across large tiled volumes.

```
HDF5 Snapshot
  └─ analyze_snapshot.py     → Delaunay tessellation, walls, filaments, clusters
      ├─ ndtopo_stats.py     → per-point topology statistics
      └─ batch_clip.py       → slab visualizations and PNG renders

batch_crop_and_clip.py       → orchestrates the above across a tiled volume
aggregate_topology_points.py → cross-crop statistics and histograms
```

## Technologies

| | |
|---|---|
| Topology extraction | [DisPerSE](https://www.iap.fr/useriap/sousbie/web/html/indexd41d.html) |
| Visualization | [ParaView](https://www.paraview.org/) / `pvpython` |
| Simulation data | [Quijote](https://quijote-simulations.readthedocs.io/en/latest/) |
| Scripting | Python |
