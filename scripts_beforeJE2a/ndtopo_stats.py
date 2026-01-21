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
    p.add_argument("--delaunay-ndnet", help="Delaunay mesh (NDnet).")
    p.add_argument("--walls-ndnet", help="Walls mesh (NDnet).")
    p.add_argument("--filaments-ndskl", help="Filaments (NDskl).")
    p.add_argument("--delaunay-vtk", help="Delaunay mesh (VTU).")
    p.add_argument("--walls-vtk", help="Walls mesh (VTU).")
    p.add_argument("--filaments-vtk", help="Filaments (VTP).")
    p.add_argument("--output-csv", required=True, help="Path to write stats CSV.")
    p.add_argument("--netconv-bin", default="netconv", help="netconv executable (for auto-conversion).")
    p.add_argument("--skelconv-bin", default="skelconv", help="skelconv executable (for auto-conversion).")
    p.add_argument("--write-vtk", action="store_true", help="Convert NDnet/NDskl to unsmoothed VTU/VTP before parsing.")
    p.add_argument("--verbose", action="store_true", help="Print progress messages.")
    p.add_argument("--delaunay-id-field", default="true_index", help="ID field for delaunay/universe (default: true_index).")
    p.add_argument(
        "--delaunay-cell-mode",
        choices=["all", "zero"],
        default="all",
        help="How to interpret delaunay cell IDs: all=int(cell) or zero=only .0 cells.",
    )
    p.add_argument("--walls-id-field", default="true_index", help="ID field for walls/universe (default: true_index).")
    p.add_argument("--filaments-id-field", default="cell", help="ID field for filaments (default: cell -> int(cell)).")
    p.add_argument(
        "--walls-cell-mode",
        choices=["all", "zero"],
        default="all",
        help="How to interpret wall cell IDs: all=int(cell) or zero=only .0 cells.",
    )
    p.add_argument(
        "--filaments-cell-mode",
        choices=["all", "zero"],
        default="all",
        help="How to interpret filament cell IDs: all=int(cell) or zero=only .0 cells.",
    )
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


def _strip_s_tag(stem: str) -> str:
    """Strip trailing smoothing tokens like .S020 or _S020 from a stem."""
    import re

    return re.sub(r"[._]S\d{3}$", "", stem)


def _sanitize_single_dot(path: Path) -> Path:
    """Ensure only the final dot remains by replacing dots in the stem with underscores."""
    clean_stem = path.stem.replace(".", "_")
    new_path = path.with_name(f"{clean_stem}{path.suffix}")
    if new_path == path:
        return path
    if new_path.exists():
        return new_path
    try:
        path.rename(new_path)
        return new_path
    except Exception:
        return path


def ensure_vtu_from_ndnet(path: Path, netconv_bin: str, role: str) -> Path:
    if path.suffix.lower() != ".ndnet":
        return path
    stem = _strip_s_tag(path.stem)
    if role == "delaunay":
        if "delaunay" not in stem:
            stem = f"{stem}_delaunay"
    elif role == "walls":
        if "manifolds" not in stem:
            stem = f"{stem}_manifolds"
    base = path.with_name(f"{stem}_S000.vtu")
    candidates = [base]
    for cand in candidates:
        if cand.exists():
            return _sanitize_single_dot(cand)
    print(f"[info] Converting NDnet to unsmoothed VTU: {path} -> {base}")
    run_cmd([netconv_bin, str(path), "-outName", base.stem, "-outDir", str(base.parent), "-to", "vtu", "-smooth", "0"])
    for cand in candidates:
        if cand.exists():
            return _sanitize_single_dot(cand)
    raise SystemExit(f"[error] Could not find converted VTU: {candidates}")


def ensure_vtp_from_ndskl(path: Path, skelconv_bin: str, role: str) -> Path:
    if path.suffix.lower() != ".ndskl":
        return path
    stem = _strip_s_tag(path.stem)
    stem = stem.replace(".", "_")
    if role == "filaments":
        tokens = stem.split("_")
        # target order: <prefix>_sX_arcs_<tag>
        if len(tokens) >= 2:
            arc_tag = tokens[-1]
            prefix_tokens = tokens[:-1]
            if "arcs" not in prefix_tokens:
                prefix_tokens.append("arcs")
            stem = "_".join(prefix_tokens + [arc_tag])
        else:
            stem = f"{stem}_arcs"
    # skelconv with -smooth 0 appends .S000; keep outName clean and sanitize after
    base = path.with_name(f"{stem}.vtp")
    candidates = [
        base,
        base.with_name(f"{stem}_S000.vtp"),
        base.with_name(f"{stem}.S000.vtp"),
    ]
    for cand in candidates:
        if cand.exists():
            return _sanitize_single_dot(cand)
    print(f"[info] Converting NDskl to unsmoothed VTP: {path} -> {base}")
    run_cmd([skelconv_bin, str(path), "-outName", stem, "-outDir", str(base.parent), "-to", "vtp", "-smooth", "0"])
    for cand in candidates:
        if cand.exists():
            return _sanitize_single_dot(cand)
    raise SystemExit(f"[error] Could not find converted VTP: {candidates}")


def read_vtk_ids(
    path: Path,
    id_field: str,
    scalar_fields: List[str],
    fallbacks: List[str],
    cell_mode: Optional[str] = None,
) -> Tuple[Dict[int, Dict[str, float]], List[int], str]:
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
    if chosen == "cell" and cell_mode == "zero":
        ids_list: List[int] = []
        for x in ids_arr:
            xf = float(x)
            if abs(xf - round(xf)) < 1e-6:
                ids_list.append(int(round(xf)))
    else:
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
        if not vals:
            out[f"{name}_sum"] = 0.0
            out[f"{name}_mean"] = 0.0
            out[f"{name}_min"] = 0.0
            out[f"{name}_q25"] = 0.0
            out[f"{name}_median"] = 0.0
            out[f"{name}_q75"] = 0.0
            out[f"{name}_max"] = 0.0
            out[f"{name}_std"] = 0.0
            continue
        vals_sorted = sorted(vals)
        n = len(vals_sorted)
        out[f"{name}_sum"] = float(sum(vals_sorted))
        out[f"{name}_mean"] = float(out[f"{name}_sum"] / n)
        out[f"{name}_min"] = float(vals_sorted[0])
        out[f"{name}_max"] = float(vals_sorted[-1])
        # simple quantiles
        def _quantile(q: float) -> float:
            if n == 1:
                return float(vals_sorted[0])
            idx = q * (n - 1)
            lo = int(idx)
            hi = min(lo + 1, n - 1)
            frac = idx - lo
            return float(vals_sorted[lo] * (1 - frac) + vals_sorted[hi] * frac)
        out[f"{name}_q25"] = _quantile(0.25)
        out[f"{name}_median"] = _quantile(0.5)
        out[f"{name}_q75"] = _quantile(0.75)
        # std
        mean_val = out[f"{name}_mean"]
        out[f"{name}_std"] = float((sum((v - mean_val) ** 2 for v in vals_sorted) / n) ** 0.5)
    return out


def main() -> None:
    args = parse_args()
    if args.verbose:
        print("[info] Inputs:")
        print(f"  delaunay: {args.delaunay_vtk or args.delaunay_ndnet}")
        print(f"  walls:    {args.walls_vtk or args.walls_ndnet}")
        print(f"  filaments:{args.filaments_vtk or args.filaments_ndskl}")
        print(
            f"  delaunay_id_field={args.delaunay_id_field}, walls_id_field={args.walls_id_field}, "
            f"filaments_id_field={args.filaments_id_field}, "
            f"delaunay_cell_mode={args.delaunay_cell_mode}, walls_cell_mode={args.walls_cell_mode}, "
            f"filaments_cell_mode={args.filaments_cell_mode}, scalars={args.scalar_fields}"
        )
    # Resolve inputs: prefer VTK if provided, otherwise require native
    delaunay_path = Path(args.delaunay_vtk) if args.delaunay_vtk else Path(args.delaunay_ndnet) if args.delaunay_ndnet else None
    walls_path = Path(args.walls_vtk) if args.walls_vtk else Path(args.walls_ndnet) if args.walls_ndnet else None
    fils_path = Path(args.filaments_vtk) if args.filaments_vtk else Path(args.filaments_ndskl) if args.filaments_ndskl else None

    if delaunay_path is None or walls_path is None or fils_path is None:
        raise SystemExit("[error] Provide either VTK or NDnet/NDskl inputs for delaunay, walls, and filaments.")

    for p in (delaunay_path, walls_path, fils_path):
        if p and not p.exists():
            raise SystemExit(f"[error] Input file not found: {p}")

    if args.write_vtk and (delaunay_path.suffix.lower() == ".ndnet" or walls_path.suffix.lower() == ".ndnet" or fils_path.suffix.lower() == ".ndskl"):
        if delaunay_path.suffix.lower() == ".ndnet":
            delaunay_path = ensure_vtu_from_ndnet(delaunay_path, args.netconv_bin, role="delaunay")
        if walls_path.suffix.lower() == ".ndnet":
            walls_path = ensure_vtu_from_ndnet(walls_path, args.netconv_bin, role="walls")
        if fils_path.suffix.lower() == ".ndskl":
            fils_path = ensure_vtp_from_ndskl(fils_path, args.skelconv_bin, role="filaments")

    id_field = args.walls_id_field
    skel_id_field = args.filaments_id_field
    scalars = args.scalar_fields
    fallback_ids = ["true_index", "index", "source_index", "cell"]

    universe, _, used_uni = read_vtk_ids(
        delaunay_path,
        args.delaunay_id_field,
        scalars,
        fallback_ids,
        cell_mode=args.delaunay_cell_mode,
    )
    _, wall_ids_list, used_w = read_vtk_ids(
        walls_path,
        id_field,
        [],
        fallback_ids,
        cell_mode=args.walls_cell_mode,
    )
    _, fil_ids_list, used_f = read_vtk_ids(
        fils_path,
        skel_id_field,
        [],
        fallback_ids,
        cell_mode=args.filaments_cell_mode,
    )
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
        "walls": wall_ids,
        "filaments": fil_ids,
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
        fieldnames.extend(
            [
                f"{name}_sum",
                f"{name}_mean",
                f"{name}_min",
                f"{name}_q25",
                f"{name}_median",
                f"{name}_q75",
                f"{name}_max",
                f"{name}_std",
            ]
        )
    with open(args.output_csv, "w", newline="", encoding="ascii") as sink:
        writer = csv.DictWriter(sink, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    print(f"[done] wrote {args.output_csv}")


if __name__ == "__main__":
    main()
