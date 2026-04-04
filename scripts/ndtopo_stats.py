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
import re
import subprocess
from pathlib import Path
from typing import Iterable


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Overlap stats between NDnet/NDskl (VTU/VTP or native).")
    p.add_argument("--delaunay-ndnet", help="Delaunay mesh (NDnet).")
    p.add_argument("--walls-ndnet", help="Walls mesh (NDnet).")
    p.add_argument("--filaments-ndskl", help="Filaments (NDskl).")
    p.add_argument("--filament-manifolds-ndnet", help="Filament manifolds (NDnet).")
    p.add_argument("--cluster-manifolds-ndnet", help="Cluster manifolds (NDnet).")
    p.add_argument("--delaunay-vtk", help="Delaunay mesh (VTU).")
    p.add_argument("--walls-vtk", help="Walls mesh (VTU).")
    p.add_argument("--filaments-vtk", help="Filaments (VTP).")
    p.add_argument("--filament-manifolds-vtk", help="Filament manifolds (VTU).")
    p.add_argument("--cluster-manifolds-vtk", help="Cluster manifolds (VTU).")
    p.add_argument("--output-csv", required=True, help="Path to write stats CSV.")
    p.add_argument(
        "--topology-scalars-csv",
        help=(
            "Optional CSV path for topology-scalar stats (scalars read from walls/filaments/"
            "filament-manifolds/clusters instead of Delaunay)."
        ),
    )
    p.add_argument(
        "--per-point-csv",
        help="Optional CSV with one row per Delaunay point (id, category, scalars).",
    )
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
    p.add_argument(
        "--filament-manifolds-id-field",
        default="true_index",
        help="ID field for filament manifolds (default: true_index).",
    )
    p.add_argument(
        "--cluster-manifolds-id-field",
        default="true_index",
        help="ID field for cluster manifolds (default: true_index).",
    )
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
        "--filament-manifolds-cell-mode",
        choices=["all", "zero"],
        default="all",
        help="How to interpret filament manifold cell IDs: all=int(cell) or zero=only .0 cells.",
    )
    p.add_argument(
        "--cluster-manifolds-cell-mode",
        choices=["all", "zero"],
        default="all",
        help="How to interpret cluster manifold cell IDs: all=int(cell) or zero=only .0 cells.",
    )
    p.add_argument(
        "--scalar-fields",
        nargs="+",
        default=["mass", "field_value", "log_field_value"],
        help="Scalar fields (from Delaunay) to sum/mean per category.",
    )
    return p.parse_args()


def run_cmd(cmd: list[str]) -> None:
    print(f"[run] {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def _strip_s_tag(stem: str) -> str:
    """Strip trailing smoothing tokens like .S020 or _S020 from a stem."""
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


def _ensure_converted(
    source: Path,
    base: Path,
    candidates: list[Path],
    cmd: list[str],
    kind: str,
) -> Path:
    """Check candidates → run converter → check again → raise on failure."""
    for cand in candidates:
        if cand.exists():
            return _sanitize_single_dot(cand)
    print(f"[info] Converting {source} -> {base}")
    run_cmd(cmd)
    for cand in candidates:
        if cand.exists():
            return _sanitize_single_dot(cand)
    raise SystemExit(f"[error] Could not find converted {kind}: {candidates}")


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
    elif role == "filament_manifolds":
        if "filament_manifolds" not in stem:
            stem = f"{stem}_filament_manifolds"
    elif role == "cluster_manifolds":
        if "cluster_manifolds" not in stem:
            stem = f"{stem}_cluster_manifolds"
    base = path.with_name(f"{stem}_S000.vtu")
    cmd = [netconv_bin, str(path), "-outName", base.stem, "-outDir", str(base.parent), "-to", "vtu", "-smooth", "0"]
    return _ensure_converted(source=path, base=base, candidates=[base], cmd=cmd, kind="VTU")


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
    cmd = [skelconv_bin, str(path), "-outName", stem, "-outDir", str(base.parent), "-to", "vtp", "-smooth", "0"]
    return _ensure_converted(source=path, base=base, candidates=candidates, cmd=cmd, kind="VTP")


def read_vtk_ids(
    path: Path,
    id_field: str,
    scalar_fields: list[str],
    fallbacks: list[str],
    cell_mode: str | None = None,
    allow_cell_fallback: bool = True,
) -> tuple[dict[int, dict[str, float]], list[int], str]:
    try:
        from vtkmodules.numpy_interface import dataset_adapter as dsa
        from vtkmodules.vtkIOXML import vtkXMLPolyDataReader, vtkXMLUnstructuredGridReader
    except Exception as exc:
        raise SystemExit(f"[error] vtkmodules not available: {exc}")

    reader = vtkXMLUnstructuredGridReader() if path.suffix.lower() == ".vtu" else vtkXMLPolyDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    dobj = dsa.WrapDataObject(reader.GetOutput())

    point_names = set(dobj.PointData.keys())
    cell_names = set(dobj.CellData.keys())
    chosen = next((c for c in [id_field] + fallbacks if c in point_names), None)
    chosen_location = "point"
    if chosen is None and allow_cell_fallback:
        chosen = next((c for c in [id_field] + fallbacks if c in cell_names), None)
        chosen_location = "cell"
    if chosen is None:
        return {}, [], id_field

    ids_arr = dobj.PointData[chosen] if chosen_location == "point" else dobj.CellData[chosen]
    if chosen == "cell" and cell_mode == "zero":
        ids_list: list[int] = []
        for x in ids_arr:
            xf = float(x)
            if abs(xf - round(xf)) < 1e-6:
                ids_list.append(int(round(xf)))
    else:
        ids_list = [int(float(x)) for x in ids_arr]

    scalars: dict[int, dict[str, float]] = {}
    if scalar_fields:
        data_block = dobj.PointData if chosen_location == "point" else dobj.CellData
        available = set(data_block.keys())
        for name in scalar_fields:
            if name not in available:
                continue
            vals = data_block[name]
            for idx, val in zip(ids_list, vals):
                scalars.setdefault(int(idx), {})[name] = float(val)
    return scalars, ids_list, chosen if chosen_location == "point" else f"{chosen} (cell)"


def list_point_arrays(path: Path) -> dict[str, list[str]]:
    try:
        from vtkmodules.numpy_interface import dataset_adapter as dsa
        from vtkmodules.vtkIOXML import vtkXMLPolyDataReader, vtkXMLUnstructuredGridReader
    except Exception:
        return {}
    reader = vtkXMLUnstructuredGridReader() if path.suffix.lower() == ".vtu" else vtkXMLPolyDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    dobj = dsa.WrapDataObject(reader.GetOutput())
    return {
        "point": list(dobj.PointData.keys()),
        "cell": list(dobj.CellData.keys()),
    }


def _quantile(vals_sorted: list[float], n: int, q: float) -> float:
    if n == 1:
        return float(vals_sorted[0])
    idx = q * (n - 1)
    lo = int(idx)
    hi = min(lo + 1, n - 1)
    frac = idx - lo
    return float(vals_sorted[lo] * (1 - frac) + vals_sorted[hi] * frac)


def aggregate(ids: Iterable[int], universe: dict[int, dict[str, float]], scalar_fields: list[str]) -> dict[str, float]:
    ids_list = list(ids)
    out: dict[str, float] = {"count": len(ids_list)}
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
        out[f"{name}_q25"] = _quantile(vals_sorted, n, 0.25)
        out[f"{name}_median"] = _quantile(vals_sorted, n, 0.5)
        out[f"{name}_q75"] = _quantile(vals_sorted, n, 0.75)
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
        if args.filament_manifolds_vtk or args.filament_manifolds_ndnet:
            print(f"  filament_manifolds:{args.filament_manifolds_vtk or args.filament_manifolds_ndnet}")
        if args.cluster_manifolds_vtk or args.cluster_manifolds_ndnet:
            print(f"  cluster_manifolds:{args.cluster_manifolds_vtk or args.cluster_manifolds_ndnet}")
        print(
            f"  delaunay_id_field={args.delaunay_id_field}, walls_id_field={args.walls_id_field}, "
            f"filaments_id_field={args.filaments_id_field}, "
            f"filament_manifolds_id_field={args.filament_manifolds_id_field}, "
            f"cluster_manifolds_id_field={args.cluster_manifolds_id_field}, "
            f"delaunay_cell_mode={args.delaunay_cell_mode}, walls_cell_mode={args.walls_cell_mode}, "
            f"filaments_cell_mode={args.filaments_cell_mode}, "
            f"filament_manifolds_cell_mode={args.filament_manifolds_cell_mode}, "
            f"cluster_manifolds_cell_mode={args.cluster_manifolds_cell_mode}, "
            f"scalars={args.scalar_fields}"
        )
    # Resolve inputs: prefer VTK if provided, otherwise require native
    delaunay_path = Path(args.delaunay_vtk) if args.delaunay_vtk else Path(args.delaunay_ndnet) if args.delaunay_ndnet else None
    walls_path = Path(args.walls_vtk) if args.walls_vtk else Path(args.walls_ndnet) if args.walls_ndnet else None
    fils_path = Path(args.filaments_vtk) if args.filaments_vtk else Path(args.filaments_ndskl) if args.filaments_ndskl else None
    filman_path = (
        Path(args.filament_manifolds_vtk)
        if args.filament_manifolds_vtk
        else Path(args.filament_manifolds_ndnet)
        if args.filament_manifolds_ndnet
        else None
    )
    cluster_path = (
        Path(args.cluster_manifolds_vtk)
        if args.cluster_manifolds_vtk
        else Path(args.cluster_manifolds_ndnet)
        if args.cluster_manifolds_ndnet
        else None
    )

    if delaunay_path is None or walls_path is None or fils_path is None:
        raise SystemExit("[error] Provide either VTK or NDnet/NDskl inputs for delaunay, walls, and filaments.")

    for p in (delaunay_path, walls_path, fils_path, filman_path, cluster_path):
        if p and not p.exists():
            raise SystemExit(f"[error] Input file not found: {p}")

    if args.write_vtk and (
        delaunay_path.suffix.lower() == ".ndnet"
        or walls_path.suffix.lower() == ".ndnet"
        or fils_path.suffix.lower() == ".ndskl"
        or (filman_path and filman_path.suffix.lower() == ".ndnet")
        or (cluster_path and cluster_path.suffix.lower() == ".ndnet")
    ):
        if delaunay_path.suffix.lower() == ".ndnet":
            delaunay_path = ensure_vtu_from_ndnet(delaunay_path, args.netconv_bin, role="delaunay")
        if walls_path.suffix.lower() == ".ndnet":
            walls_path = ensure_vtu_from_ndnet(walls_path, args.netconv_bin, role="walls")
        if fils_path.suffix.lower() == ".ndskl":
            fils_path = ensure_vtp_from_ndskl(fils_path, args.skelconv_bin, role="filaments")
        if filman_path and filman_path.suffix.lower() == ".ndnet":
            filman_path = ensure_vtu_from_ndnet(filman_path, args.netconv_bin, role="filament_manifolds")
        if cluster_path and cluster_path.suffix.lower() == ".ndnet":
            cluster_path = ensure_vtu_from_ndnet(cluster_path, args.netconv_bin, role="cluster_manifolds")

    id_field = args.walls_id_field
    skel_id_field = args.filaments_id_field
    scalars = args.scalar_fields
    fallback_ids = ["true_index", "index", "source_index", "cell"]

    need_topo_scalars = bool(args.topology_scalars_csv)
    universe, uni_ids_list, delaunay_id_used = read_vtk_ids(
        delaunay_path,
        args.delaunay_id_field,
        scalars,
        fallback_ids,
        cell_mode=args.delaunay_cell_mode,
        allow_cell_fallback=False,
    )
    wall_scalars, wall_ids_list, walls_id_used = read_vtk_ids(
        walls_path,
        id_field,
        scalars if need_topo_scalars else [],
        fallback_ids,
        cell_mode=args.walls_cell_mode,
    )
    fil_scalars, fil_ids_list, fils_id_used = read_vtk_ids(
        fils_path,
        skel_id_field,
        scalars if need_topo_scalars else [],
        fallback_ids,
        cell_mode=args.filaments_cell_mode,
    )
    filman_ids_list: list[int] = []
    filman_scalars: dict[int, dict[str, float]] = {}
    filman_id_used = ""
    if filman_path is not None:
        filman_scalars, filman_ids_list, filman_id_used = read_vtk_ids(
            filman_path,
            args.filament_manifolds_id_field,
            scalars if need_topo_scalars else [],
            fallback_ids,
            cell_mode=args.filament_manifolds_cell_mode,
        )
    cluster_ids_list: list[int] = []
    cluster_scalars: dict[int, dict[str, float]] = {}
    cluster_id_used = ""
    if cluster_path is not None:
        cluster_scalars, cluster_ids_list, cluster_id_used = read_vtk_ids(
            cluster_path,
            args.cluster_manifolds_id_field,
            scalars if need_topo_scalars else [],
            fallback_ids,
            cell_mode=args.cluster_manifolds_cell_mode,
        )
    wall_ids = set(wall_ids_list)
    fil_ids = set(fil_ids_list)
    filman_ids = set(filman_ids_list)
    cluster_ids = set(cluster_ids_list)

    if args.verbose:
        extra = f" filman:{len(filman_ids)} ({filman_id_used})" if filman_path is not None else ""
        if cluster_path is not None:
            extra += f" cluster:{len(cluster_ids)} ({cluster_id_used})"
        print(
            f"[info] Loaded ids -> universe:{len(universe)} ({delaunay_id_used}) walls:{len(wall_ids)} ({walls_id_used}) "
            f"filaments:{len(fil_ids)} ({fils_id_used}){extra}"
        )
    if len(universe) == 0 or len(wall_ids) == 0 or len(fil_ids) == 0:
        raise SystemExit("[error] Unable to load IDs from provided inputs. Check id fields and file formats.")

    uni_ids = set(uni_ids_list or universe.keys())

    if cluster_path is not None and not cluster_ids:
        cluster_id_used_plain = cluster_id_used.replace(" (cell)", "")
        fallback_candidates = [name for name in ("cell", "source_index", "index", "true_index") if name != cluster_id_used_plain]
        for candidate in fallback_candidates:
            alt_scalars, alt_ids_list, alt_used = read_vtk_ids(
                cluster_path,
                candidate,
                scalars if need_topo_scalars else [],
                fallback_ids,
                cell_mode=args.cluster_manifolds_cell_mode,
            )
            if alt_ids_list:
                print(
                    f"[warn] Cluster IDs empty with '{cluster_id_used}', falling back to '{alt_used}' for {cluster_path}."
                )
                cluster_ids = set(alt_ids_list)
                cluster_id_used = alt_used
                if need_topo_scalars:
                    cluster_scalars = alt_scalars
                break
        if not cluster_ids:
            arrays = list_point_arrays(cluster_path)
            if arrays:
                print(
                    f"[warn] Cluster IDs still empty for {cluster_path}. "
                    f"Available arrays: point={arrays.get('point', [])}, cell={arrays.get('cell', [])}. "
                    "Try --cluster-manifolds-id-field <name> (often 'cell' or 'source_index')."
                )
            else:
                print(
                    f"[warn] Cluster IDs still empty for {cluster_path}. "
                    "Try --cluster-manifolds-id-field <name> (often 'cell' or 'source_index')."
                )
    shared_walls_filaments = wall_ids & fil_ids
    walls_not_filaments = wall_ids - fil_ids
    filaments_not_walls = fil_ids - wall_ids
    unassigned = uni_ids - (wall_ids | fil_ids | filman_ids | cluster_ids)

    category_list: list[tuple[str, set[int]]] = [
        ("walls", wall_ids),
        ("filaments", fil_ids),
        ("walls_not_filaments", walls_not_filaments),
        ("filaments_not_walls", filaments_not_walls),
        ("shared_walls_filaments", shared_walls_filaments),
    ]

    if filman_path is not None:
        shared_walls_filament_manifolds = wall_ids & filman_ids
        shared_filaments_filament_manifolds = fil_ids & filman_ids
        walls_not_filament_manifolds = wall_ids - filman_ids
        filament_manifolds_not_walls = filman_ids - wall_ids
        category_list.extend(
            [
                ("filament_manifolds", filman_ids),
                ("walls_not_filament_manifolds", walls_not_filament_manifolds),
                ("filament_manifolds_not_walls", filament_manifolds_not_walls),
                ("shared_walls_filament_manifolds", shared_walls_filament_manifolds),
                ("shared_filaments_filament_manifolds", shared_filaments_filament_manifolds),
            ]
        )

    if cluster_path is not None:
        clusters = cluster_ids
        clusters_not_filaments = clusters - fil_ids
        clusters_not_walls = clusters - wall_ids
        filaments_not_clusters = fil_ids - clusters
        walls_not_clusters = wall_ids - clusters
        shared_walls_clusters = wall_ids & clusters
        shared_filaments_clusters = fil_ids & clusters
        shared_walls_filaments_clusters = wall_ids & fil_ids & clusters
        category_list.extend(
            [
                ("clusters", clusters),
                ("clusters_not_filaments", clusters_not_filaments),
                ("clusters_not_walls", clusters_not_walls),
                ("filaments_not_clusters", filaments_not_clusters),
                ("walls_not_clusters", walls_not_clusters),
                ("shared_walls_clusters", shared_walls_clusters),
                ("shared_filaments_clusters", shared_filaments_clusters),
                ("shared_walls_filaments_clusters", shared_walls_filaments_clusters),
            ]
        )
        if filman_path is not None:
            clusters_not_filament_manifolds = clusters - filman_ids
            filament_manifolds_not_clusters = filman_ids - clusters
            shared_filament_manifolds_clusters = filman_ids & clusters
            shared_walls_filament_manifolds_clusters = wall_ids & filman_ids & clusters
            shared_filaments_filament_manifolds_clusters = fil_ids & filman_ids & clusters
            shared_walls_filaments_filament_manifolds_clusters = (
                wall_ids & fil_ids & filman_ids & clusters
            )
            category_list.extend(
                [
                    ("clusters_not_filament_manifolds", clusters_not_filament_manifolds),
                    ("filament_manifolds_not_clusters", filament_manifolds_not_clusters),
                    ("shared_filament_manifolds_clusters", shared_filament_manifolds_clusters),
                    (
                        "shared_walls_filament_manifolds_clusters",
                        shared_walls_filament_manifolds_clusters,
                    ),
                    (
                        "shared_filaments_filament_manifolds_clusters",
                        shared_filaments_filament_manifolds_clusters,
                    ),
                    (
                        "shared_walls_filaments_filament_manifolds_clusters",
                        shared_walls_filaments_filament_manifolds_clusters,
                    ),
                ]
            )

    category_list.append(("unassigned", unassigned))

    rows: list[dict[str, object]] = []
    for label, ids in category_list:
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

    if args.topology_scalars_csv:
        sources: list[tuple[str, dict[int, dict[str, float]]]] = []
        if wall_scalars:
            sources.append(("walls", wall_scalars))
        if fil_scalars:
            sources.append(("filaments", fil_scalars))
        if filman_path is not None and filman_scalars:
            sources.append(("filament_manifolds", filman_scalars))
        if cluster_path is not None and cluster_scalars:
            sources.append(("clusters", cluster_scalars))
        if not sources:
            print("[warn] topology-scalars requested but no topology scalar fields were found.")
        else:
            topo_fieldnames = ["category", "scalar_source", "count"]
            for name in scalars:
                topo_fieldnames.extend(
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
            with open(args.topology_scalars_csv, "w", newline="", encoding="ascii") as sink:
                writer = csv.DictWriter(sink, fieldnames=topo_fieldnames)
                writer.writeheader()
                for label, ids in category_list:
                    for source_name, source_scalars in sources:
                        source_ids = set(source_scalars.keys())
                        ids_for_source = source_ids & ids
                        stats = aggregate(ids_for_source, source_scalars, scalars)
                        stats["category"] = label
                        stats["scalar_source"] = source_name
                        writer.writerow(stats)
            print(f"[done] wrote {args.topology_scalars_csv}")

    if args.per_point_csv:
        per_point_fields = [
            "delaunay_id",
            "is_wall",
            "is_filament",
            "is_filament_manifold",
            "is_cluster",
        ] + scalars
        with open(args.per_point_csv, "w", newline="", encoding="ascii") as sink:
            writer = csv.DictWriter(sink, fieldnames=per_point_fields)
            writer.writeheader()
            for idx in sorted(uni_ids):
                in_wall = idx in wall_ids
                in_fil = idx in fil_ids
                in_filman = filman_path is not None and idx in filman_ids
                in_cluster = cluster_path is not None and idx in cluster_ids
                row = {
                    "delaunay_id": idx,
                    "is_wall": int(in_wall),
                    "is_filament": int(in_fil),
                    "is_filament_manifold": int(in_filman),
                    "is_cluster": int(in_cluster),
                }
                vals = universe.get(idx, {})
                for name in scalars:
                    row[name] = vals.get(name, 0.0)
                writer.writerow(row)
        print(f"[done] wrote {args.per_point_csv}")


if __name__ == "__main__":
    main()
