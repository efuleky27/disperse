#!/usr/bin/env python3
"""Assign each (sampled) particle to a basin by nearest minimum critical point.

This is a lightweight, post-processing helper: given a particle cloud (VTU)
and a critical-points dump (crits_ascii) from DisPerSE, it attaches a
`basin_id` array to every particle based on the nearest minimum. This is a
proximity-based proxy for basin ownership; it does not traverse the Morse–Smale
complex itself. For exact Morse–Smale cell ownership, a dedicated NDnet
segmentation step would be required.

Usage
-----
1) Export your particle cloud to VTU (matching the decimation/crop used in
   DisPerSE), e.g.:
   python scripts/export_snapshot_vtu.py --input data/snap_010.hdf5 \
       --parttype PartType1 --output outputs/snap_010_points.vtu \
       --target-count 2000000

2) Dump critical points to ASCII (minima included) using skelconv:
   skelconv your_skeleton.NDskl -to crits_ascii -outDir outputs -outName crits
   This produces e.g. outputs/crits.crits

3) Run the basin labeler:
   python scripts/label_basin_assignment.py \
       --points-vtu outputs/snap_010_points.vtu \
       --crits-ascii outputs/crits.crits \
       --output outputs/snap_010_points_labeled.vtu

The output VTU will carry two new point-data arrays:
  - basin_id (int): index of the nearest minimum
  - basin_distance (float): distance to that minimum (same units as the VTU)
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple

import numpy as np
import vtk
from scipy.spatial import cKDTree
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy


@dataclass
class CriticalPoint:
    idx: int
    kind: str
    pos: np.ndarray


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Assign each particle to the nearest minimum critical point (basin proxy)."
    )
    parser.add_argument("--points-vtu", required=True, help="VTU file containing particle points.")
    parser.add_argument(
        "--crits-ascii",
        required=True,
        help="ASCII critical points from skelconv -to crits_ascii (expects minima present).",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output VTU with basin_id and basin_distance point-data arrays.",
    )
    parser.add_argument(
        "--min-kind",
        default="0",
        help="Value in the crits file representing minima (default: '0'). Accepts string tokens too (e.g., MIN).",
    )
    return parser.parse_args()


def load_critical_points(path: Path, min_kind: str) -> List[CriticalPoint]:
    """Parse crits_ascii (id kind x y z ...) and keep minima."""
    minima: List[CriticalPoint] = []
    raw = Path(path).read_text().strip().splitlines()
    for line in raw:
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 5:
            continue
        idx = int(parts[0])
        kind = parts[1]
        if kind != min_kind:
            continue
        try:
            x, y, z = map(float, parts[2:5])
        except Exception:
            continue
        minima.append(CriticalPoint(idx=idx, kind=kind, pos=np.array([x, y, z], dtype=np.float64)))
    if not minima:
        raise SystemExit(f"No minima with kind '{min_kind}' found in {path}")
    return minima


def load_points(path: Path) -> Tuple[np.ndarray, vtk.vtkUnstructuredGrid]:
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    data = reader.GetOutput()
    pts_vtk = data.GetPoints().GetData()
    pts = vtk_to_numpy(pts_vtk).astype(np.float64, copy=False)
    return pts, data


def attach_array(grid: vtk.vtkUnstructuredGrid, name: str, arr: np.ndarray) -> None:
    vtk_arr = numpy_to_vtk(arr, deep=True)
    vtk_arr.SetName(name)
    grid.GetPointData().AddArray(vtk_arr)


def main() -> None:
    args = parse_args()
    points_path = Path(args.points_vtu)
    crits_path = Path(args.crits_ascii)
    output_path = Path(args.output)

    minima = load_critical_points(crits_path, args.min_kind)
    min_positions = np.vstack([c.pos for c in minima])
    min_ids = np.array([c.idx for c in minima], dtype=np.int64)

    pts, grid = load_points(points_path)
    tree = cKDTree(min_positions)
    distances, nn_idx = tree.query(pts)
    basin_ids = min_ids[nn_idx]

    attach_array(grid, "basin_id", basin_ids.astype(np.int64))
    attach_array(grid, "basin_distance", distances.astype(np.float64))

    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(str(output_path))
    writer.SetInputData(grid)
    writer.SetDataModeToBinary()
    writer.Write()
    print(f"[done] Labeled {pts.shape[0]} points -> {output_path}")


if __name__ == "__main__":
    main()
