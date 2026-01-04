#!/usr/bin/env python3
"""
Compute overlap/unique/unassigned stats between a Delaunay mesh (VTU),
manifold mesh (walls, VTU), and skeleton (filaments, VTP). IDs are matched
between files; scalars are always taken from the Delaunay mesh for aggregation.

If you pass native NDnet/NDskl, use --write-vtk to auto-convert them to
unsmoothed VTU/VTP (netconv/skelconv with -smooth 0) before processing.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Overlap stats between NDnet/NDskl (VTU/VTP or native).")
    p.add_argument("--delaunay-ndnet", required=True, help="Delaunay mesh (VTU/NDnet).")
    p.add_argument("--walls-ndnet", required=True, help="Walls mesh (VTU/NDnet).")
    p.add_argument("--filaments-ndskl", required=True, help="Filaments (VTP/NDskl).")
    p.add_argument("--output-csv", required=True, help="Path to write stats CSV.")
    p.add_argument("--netconv-bin", default="netconv", help="netconv executable (for auto-conversion).")
    p.add_argument("--skelconv-bin", default="skelconv", help="skelconv executable (for auto-conversion).")
    p.add_argument("--write-vtk", action="store_true", help="Convert NDnet/NDskl to unsmoothed VTU/VTP before parsing.")
    p.add_argument("--verbose", action="store_true", help="Print progress messages.")
    p.add_argument("--walls-id-field", default="true_index", help="ID field for walls/universe (default: true_index).")
    p.add_argument("--filaments-id-field", default="cell", help="ID field for filaments (default: cell -> int(cell)).")
    p.add_argument(
        "--scalar-fields",
        nargs="+",
        default=["mass", "field_value", "log_field_value"],
        help="Scalar fields (from Delaunay) to sum/mean per category.",
    )
    return p.parse_args()


def run_cmd(cmd: List[str]) -> None:
    print(f"[run] {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def ensure_vtu_from_ndnet(path: Path, netconv_bin: str) -> Path:
    if path.suffix.lower() != ".ndnet":
        return path
    base = path.with_name(f"{path.stem}_nosmooth.vtu")
    candidates = [base, base.with_name(f"{base.stem}.S000{base.suffix}"), base.with_name(f"{base.stem}_S000{base.suffix}")]
    for cand in candidates:
        if cand.exists():
            return cand
    print(f"[info] Converting NDnet to unsmoothed VTU: {path} -> {base}")
    run_cmd([netconv_bin, str(path), "-outName", base.stem, "-outDir", str(base.parent), "-to", "vtu", "-smooth", "0"])
    for cand in candidates:
        if cand.exists():
            return cand
    raise SystemExit(f"[error] Could not find converted VTU: {candidates}")


def ensure_vtp_from_ndskl(path: Path, skelconv_bin: str) -> Path:
    if path.suffix.lower() != ".ndskl":
        return path
    base = path.with_name(f"{path.stem}_nosmooth.vtp")
    candidates = [base, base.with_name(f"{base.stem}.S000{base.suffix}"), base.with_name(f"{base.stem}_S000{base.suffix}")]
    for cand in candidates:
        if cand.exists():
            return cand
    print(f"[info] Converting NDskl to unsmoothed VTP: {path} -> {base}")
    run_cmd([skelconv_bin, str(path), "-outName", base.stem, "-outDir", str(base.parent), "-to", "vtp", "-smooth", "0"])
    for cand in candidates:
        if cand.exists():
            return cand
    raise SystemExit(f"[error] Could not find converted VTP: {candidates}")


def read_vtk_ids(path: Path, id_field: str, scalar_fields: List[str], fallbacks: List[str]) -> Tuple[Dict[int, Dict[str, float]], List[int], str]:
    try:
        from vtkmodules.numpy_interface import dataset_adapter as dsa
        from vtkmodules.vtkIOXML import vtkXMLPolyDataReader, vtkXMLUnstructuredGridReader
    except Exception as exc:
        raise SystemExit(f"[error] vtkmodules not available: {exc}")

    reader = vtkXMLUnstructuredGridReader() if path.suffix.lower() == ".vtu" else vtkXMLPolyDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    dobj = dsa.WrapDataObject(reader.GetOutput())

    chosen = next((c for c in [id_field] + fallbacks if c in dobj.PointData.keys()), None)
    if chosen is None:
        return {}, [], id_field

    ids_arr = dobj.PointData[chosen]
    ids_list = [int(float(x)) for x in ids_arr]

    scalars: Dict[int, Dict[str, float]] = {}
    for name in scalar_fields:
        if name not in dobj.PointData.keys():
            continue
        vals = dobj.PointData[name]
        for idx, val in zip(ids_list, vals):
            scalars.setdefault(int(idx), {})[name] = float(val)
    return scalars, ids_list, chosen


def aggregate(ids: Iterable[int], universe: Dict[int, Dict[str, float]], scalar_fields: List[str]) -> Dict[str, float]:
    ids_list = list(ids)
    out: Dict[str, float] = {"count": len(ids_list)}
    for name in scalar_fields:
        vals = [universe[i][name] for i in ids_list if i in universe and name in universe[i]]
        out[f"{name}_sum"] = float(sum(vals)) if vals else 0.0
        out[f"{name}_mean"] = float(sum(vals) / len(vals)) if vals else 0.0
    return out


def main() -> None:
    args = parse_args()
    if args.verbose:
        print("[info] Inputs:")
        print(f"  delaunay: {args.delaunay_ndnet}")
        print(f"  walls:    {args.walls_ndnet}")
        print(f"  filaments:{args.filaments_ndskl}")
        print(f"  walls_id_field={args.walls_id_field}, filaments_id_field={args.filaments_id_field}, scalars={args.scalar_fields}")
    for p in (args.delaunay_ndnet, args.walls_ndnet, args.filaments_ndskl):
        if p and not Path(p).exists():
            raise SystemExit(f"[error] Input file not found: {p}")

    delaunay_path = Path(args.delaunay_ndnet)
    walls_path = Path(args.walls_ndnet)
    fils_path = Path(args.filaments_ndskl)
    if args.write_vtk:
        delaunay_path = ensure_vtu_from_ndnet(delaunay_path, args.netconv_bin)
        walls_path = ensure_vtu_from_ndnet(walls_path, args.netconv_bin)
        fils_path = ensure_vtp_from_ndskl(fils_path, args.skelconv_bin)

    id_field = args.walls_id_field
    skel_id_field = args.filaments_id_field
    scalars = args.scalar_fields
    fallback_ids = ["true_index", "index", "source_index", "cell"]

    universe, _, used_uni = read_vtk_ids(delaunay_path, id_field, scalars, fallback_ids)
    _, wall_ids_list, used_w = read_vtk_ids(walls_path, id_field, [], fallback_ids)
    _, fil_ids_list, used_f = read_vtk_ids(fils_path, skel_id_field, [], fallback_ids)
    wall_ids = set(wall_ids_list)
    fil_ids = set(fil_ids_list)

    if args.verbose:
        print(f"[info] Loaded ids -> universe:{len(universe)} ({used_uni}) walls:{len(wall_ids)} ({used_w}) filaments:{len(fil_ids)} ({used_f})")
    if len(universe) == 0 or len(wall_ids) == 0 or len(fil_ids) == 0:
        raise SystemExit("[error] Unable to load IDs from provided inputs. Check id fields and file formats.")

    uni_ids = set(universe.keys())
    shared = wall_ids & fil_ids
    wall_only = wall_ids - fil_ids
    fil_only = fil_ids - wall_ids
    unassigned = uni_ids - (wall_ids | fil_ids)

    categories = {
        "shared": shared,
        "wall_only": wall_only,
        "filament_only": fil_only,
        "unassigned": unassigned,
    }

    rows: List[Dict[str, object]] = []
    for label, ids in categories.items():
        stats = aggregate(ids, universe, scalars)
        stats["category"] = label
        rows.append(stats)

    fieldnames = ["category", "count"]
    for name in scalars:
        fieldnames.extend([f"{name}_sum", f"{name}_mean"])
    with open(args.output_csv, "w", newline="", encoding="ascii") as sink:
        writer = csv.DictWriter(sink, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    print(f"[done] wrote {args.output_csv}")


if __name__ == "__main__":
    main()
