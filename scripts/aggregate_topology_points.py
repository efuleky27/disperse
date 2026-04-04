#!/usr/bin/env python3
"""
Aggregate per-point topology CSVs across crops into summary stats and histograms.

Reads *_topology_points.csv files (from ndtopo_stats.py --per-point-csv) and
computes the same category stats as topology_stats.csv, plus histogram outputs.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import math
from pathlib import Path
from collections.abc import Iterable


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Aggregate topology_points.csv files across crops.")
    p.add_argument("--root", default="outputs", help="Root directory to search for topology_points.csv files.")
    p.add_argument("--glob", default="**/*_topology_points.csv", help="Glob pattern under --root.")
    p.add_argument("--inputs", nargs="+", help="Explicit list of topology_points.csv files.")
    p.add_argument("--output-dir", required=True, help="Directory for combined outputs.")
    p.add_argument("--output-prefix", default="topology_combined", help="Prefix for output files.")
    p.add_argument("--bins", type=int, default=100, help="Histogram bins per scalar.")
    p.add_argument(
        "--hist-scalars",
        nargs="+",
        help="Scalars to histogram (default: field_value log_field_value). Use 'all' for every scalar.",
    )
    p.add_argument(
        "--hist-bin-mode",
        choices=["per-category", "global"],
        default="per-category",
        help="Histogram binning mode (default: per-category). 'global' uses shared bins per scalar.",
    )
    p.add_argument(
        "--engine",
        choices=["polars", "python", "stream"],
        default="polars",
        help="Aggregation engine: polars (fast), stream (approx), python (exact, slower).",
    )
    p.add_argument(
        "--polars-chunks",
        type=int,
        default=1,
        help="Split inputs into N chunks when using polars to reduce memory (default: 1).",
    )
    p.add_argument(
        "--hist-percentile-range",
        type=float,
        nargs=2,
        metavar=("PLOW", "PHIGH"),
        help="Percentile range for histogram bounds (e.g., 1 99).",
    )
    p.add_argument(
        "--plot-percentile-range",
        type=float,
        nargs=2,
        metavar=("PLOW", "PHIGH"),
        help="Percentile range to trim per-category data before violin/box plots (e.g., 1 99).",
    )
    p.add_argument(
        "--violin-scalar",
        default="log_field_value",
        help="Scalar to display in the violin plot (default: log_field_value).",
    )
    p.add_argument(
        "--box-scalar",
        default="field_value",
        help="Scalar to display in the box plot (default: field_value).",
    )
    p.add_argument("--no-plots", action="store_true", help="Skip PNG histogram plots.")
    p.add_argument(
        "--plots-from-raw",
        action="store_true",
        help="Build violin/box plots from raw per-point values (slower).",
    )
    p.add_argument(
        "--plot-sample-size",
        type=int,
        default=1_000_000_000,
        help="Max samples per category when building plots from histograms.",
    )
    p.add_argument(
        "--plot-fontscale",
        type=float,
        default=1.0,
        help="Scale factor for plot title/label/table fonts (default: 1.0).",
    )
    p.add_argument(
        "--plot-dpi",
        type=int,
        default=200,
        help="DPI for saved PNG plots (default: 200).",
    )
    p.add_argument(
        "--log10-field-value",
        action="store_true",
        help="Convert field_value to log10 before all aggregation and plotting.",
    )
    return p.parse_args()


def _quantile(vals_sorted: list[float], q: float) -> float:
    n = len(vals_sorted)
    if n == 1:
        return float(vals_sorted[0])
    idx = q * (n - 1)
    lo = int(idx)
    hi = min(lo + 1, n - 1)
    frac = idx - lo
    return float(vals_sorted[lo] * (1 - frac) + vals_sorted[hi] * frac)


def _maybe_log10_field_value(val: float, enabled: bool) -> float | None:
    if not enabled:
        return val
    if val <= 0:
        return None
    return math.log10(val)




def _select_hist_scalars(all_scalars: list[str], args: argparse.Namespace) -> list[str]:
    if args.hist_scalars:
        if len(args.hist_scalars) == 1 and args.hist_scalars[0].lower() == "all":
            scalars = list(all_scalars)
        else:
            scalars = [name for name in args.hist_scalars if name in all_scalars]
    else:
        defaults = ["field_value", "log_field_value"]
        scalars = [name for name in defaults if name in all_scalars]
    if args.log10_field_value and "log10_field_value" in all_scalars:
        scalars = [name for name in scalars if name != "field_value"]
        scalars.append("log10_field_value")
    return scalars


def _font_sizes(scale: float) -> dict[str, float]:
    base = {"title": 12.0, "label": 11.0, "table": 10.0, "tick": 10.0}
    return {key: val * scale for key, val in base.items()}


def _scalar_label(name: str) -> str:
    if name == "field_value":
        return "Density"
    if name == "log_field_value":
        return "$Ln$-Density"
    if name == "log10_field_value":
        return "$Log_{10}$-Density"
    return name


def _scalar_title(name: str) -> str:
    if name == "field_value":
        return "Density"
    if name == "log_field_value":
        return "Ln-Density"
    if name == "log10_field_value":
        return "Log10-Density"
    return name


def _plot_scalar_name(name: str, log10_field_value: bool) -> str:
    if log10_field_value and name == "field_value":
        return "log10_field_value"
    return name


def _apply_text_outline(text_obj, linewidth: float = 1.0) -> None:
    try:
        import matplotlib.patheffects as pe
    except Exception:
        return
    text_obj.set_path_effects([pe.withStroke(linewidth=linewidth, foreground="black")])


def _center_pad(text: str, width: int) -> str:
    figure_space = "\u2007"
    figure_dash = "\u2012"
    if text.startswith("-"):
        text = figure_dash + text[1:]
    if len(text) >= width:
        return text
    pad = width - len(text)
    left = pad // 2
    right = pad - left
    return f"{figure_space * left}{text}{figure_space * right}"


def _fmt_table_num(val: float, precision: int = 3, width: int = 9) -> str:
    # Center padded string so decimal appears near the middle of the cell.
    return _center_pad(f"{val:.{precision}f}", width)


def _fmt_table_sci(val: float, precision: int = 2, width: int = 9) -> str:
    return _center_pad(f"{val:.{precision}e}", width)


def _split_chunks(items: list[Path], chunks: int) -> list[list[Path]]:
    if chunks <= 1:
        return [items]
    buckets: list[list[Path]] = [[] for _ in range(chunks)]
    for idx, item in enumerate(items):
        buckets[idx % chunks].append(item)
    return [b for b in buckets if b]



def _composition_lists() -> tuple[list[str], list[str], list[str]]:
    keys: list[str] = []
    labels: list[str] = []
    colors: list[str] = []
    for key, label, color in COMPOSITION_COMPONENTS:
        keys.append(key)
        labels.append(label)
        colors.append(color)
    return keys, labels, colors


def _composition_color_map() -> dict[int, str]:
    return {idx: color for idx, (_, _, color) in enumerate(COMPOSITION_COMPONENTS)}


def _composition_counts_from_rows(stats_rows: list[dict[str, object]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in stats_rows:
        cat = row.get("category")
        if not cat:
            continue
        try:
            counts[str(cat)] = int(row.get("count", 0) or 0)
        except (TypeError, ValueError):
            counts[str(cat)] = 0
    return counts


def _composition_sums_from_rows(
    stats_rows: list[dict[str, object]], scalar_name: str
) -> dict[str, float]:
    sums: dict[str, float] = {}
    key = f"{scalar_name}_sum"
    for row in stats_rows:
        cat = row.get("category")
        if not cat:
            continue
        try:
            sums[str(cat)] = float(row.get(key, 0.0) or 0.0)
        except (TypeError, ValueError):
            sums[str(cat)] = 0.0
    return sums


def _seed_category_rows(
    rows: list[dict[str, object]], scalars: list[str], categories: list[str]
) -> list[dict[str, object]]:
    existing = {row.get("category") for row in rows}
    seeded = list(rows)
    for cat in categories:
        if cat in existing:
            continue
        row: dict[str, object] = {"category": cat, "count": 0}
        for name in scalars:
            for key in ("sum", "mean", "min", "q25", "median", "q75", "max", "std"):
                row[f"{name}_{key}"] = 0.0
        seeded.append(row)
    return seeded


def _write_transform_stats(
    out_dir: Path,
    output_prefix: str,
    stats_rows: list[dict[str, object]],
    scalars: list[str],
) -> None:
    out_path = out_dir / f"{output_prefix}_topology_stats_transforms.csv"
    fieldnames = [
        "category",
        "transform",
        "scalar",
        "count",
        "sum",
        "mean",
        "min",
        "q25",
        "median",
        "q75",
        "max",
        "std",
    ]
    with open(out_path, "w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in stats_rows:
            category = row.get("category", "")
            count = row.get("count", 0)
            for name in scalars:
                writer.writerow(
                    {
                        "category": category,
                        "transform": _scalar_title(name),
                        "scalar": name,
                        "count": count,
                        "sum": row.get(f"{name}_sum", 0.0),
                        "mean": row.get(f"{name}_mean", 0.0),
                        "min": row.get(f"{name}_min", 0.0),
                        "q25": row.get(f"{name}_q25", 0.0),
                        "median": row.get(f"{name}_median", 0.0),
                        "q75": row.get(f"{name}_q75", 0.0),
                        "max": row.get(f"{name}_max", 0.0),
                        "std": row.get(f"{name}_std", 0.0),
                    }
                )
    print(f"[done] wrote {out_path}")


def _infer_schema_from_headers(paths: list[Path]) -> tuple[list[str], bool, bool]:
    scalars: list[str] = []
    include_filman = False
    include_cluster = False
    skip_cols = {
        "delaunay_id",
        "is_wall",
        "is_filament",
        "is_filament_manifold",
        "is_cluster",
        "is_cluster_manifold",
    }
    scalar_set = set()
    for path in paths:
        with open(path, newline="", encoding="ascii") as handle:
            reader = csv.DictReader(handle)
            fields = reader.fieldnames or []
            if "is_filament_manifold" in fields:
                include_filman = True
            if "is_cluster" in fields or "is_cluster_manifold" in fields:
                include_cluster = True
            for name in fields:
                if name and name not in skip_cols:
                    scalar_set.add(name)
    scalars = sorted(scalar_set)
    if "field_value" in scalars and "log10_field_value" not in scalars:
        scalars.append("log10_field_value")
    return scalars, include_filman, include_cluster


def _build_fieldnames(scalars: list[str]) -> list[str]:
    fieldnames = ["category", "count"]
    for name in scalars:
        fieldnames.extend([
            f"{name}_sum", f"{name}_mean", f"{name}_min",
            f"{name}_q25", f"{name}_median", f"{name}_q75",
            f"{name}_max", f"{name}_std",
        ])
    return fieldnames


def _set_transparent(fig, ax) -> None:
    try:
        fig.patch.set_alpha(0.0)
    except Exception as exc:
        print(f"[warn] Could not set figure background alpha: {exc}")
    try:
        ax.patch.set_alpha(0.0)
    except Exception as exc:
        print(f"[warn] Could not set axes background alpha: {exc}")


def aggregate_values(vals: list[float]) -> dict[str, float]:
    if not vals:
        return {
            "count": 0,
            "sum": 0.0,
            "mean": 0.0,
            "min": 0.0,
            "q25": 0.0,
            "median": 0.0,
            "q75": 0.0,
            "max": 0.0,
            "std": 0.0,
            "skew": 0.0,
        }
    vals_sorted = sorted(vals)
    n = len(vals_sorted)
    total = float(sum(vals_sorted))
    mean_val = total / n
    std_val = float((sum((v - mean_val) ** 2 for v in vals_sorted) / n) ** 0.5)
    if std_val > 0:
        median_val = _quantile(vals_sorted, 0.5)
        skew_val = float((mean_val - median_val) / std_val)
    else:
        skew_val = 0.0
    return {
        "count": n,
        "sum": total,
        "mean": mean_val,
        "min": float(vals_sorted[0]),
        "q25": _quantile(vals_sorted, 0.25),
        "median": _quantile(vals_sorted, 0.5),
        "q75": _quantile(vals_sorted, 0.75),
        "max": float(vals_sorted[-1]),
        "std": std_val,
        "skew": skew_val,
    }


def aggregate_values_np(arr) -> dict[str, float]:
    try:
        import numpy as np
    except Exception:
        return aggregate_values(list(arr))
    if arr is None:
        return {
            "count": 0,
            "sum": 0.0,
            "mean": 0.0,
            "min": 0.0,
            "q25": 0.0,
            "median": 0.0,
            "q75": 0.0,
            "max": 0.0,
            "std": 0.0,
            "skew": 0.0,
        }
    arr_np = np.asarray(arr, dtype=float)
    if arr_np.size == 0:
        return {
            "count": 0,
            "sum": 0.0,
            "mean": 0.0,
            "min": 0.0,
            "q25": 0.0,
            "median": 0.0,
            "q75": 0.0,
            "max": 0.0,
            "std": 0.0,
            "skew": 0.0,
        }
    count = int(arr_np.size)
    total = float(arr_np.sum())
    mean_val = total / count
    q25, median, q75 = np.quantile(arr_np, [0.25, 0.5, 0.75])
    std_val = float(arr_np.std())
    if std_val > 0:
        median_val = float(np.median(arr_np))
        skew_val = float((mean_val - median_val) / std_val)
    else:
        skew_val = 0.0
    return {
        "count": count,
        "sum": total,
        "mean": mean_val,
        "min": float(arr_np.min()),
        "q25": float(q25),
        "median": float(median),
        "q75": float(q75),
        "max": float(arr_np.max()),
        "std": std_val,
        "skew": skew_val,
    }


def trim_by_percentile(values, plow: float | None, phigh: float | None):
    if values is None:
        return values
    if plow is None or phigh is None:
        return values
    try:
        import numpy as np
    except Exception:
        return values
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return arr
    qlo = plow / 100.0
    qhi = phigh / 100.0
    lo, hi = np.quantile(arr, [qlo, qhi])
    return arr[(arr >= lo) & (arr <= hi)]


def _normalize_edges(edges: list[float]) -> list[float]:
    if not edges:
        return edges
    out = [float(edges[0])]
    eps = 1e-12
    for i in range(1, len(edges)):
        val = float(edges[i])
        if val <= out[-1]:
            val = out[-1] + eps
        out.append(val)
    return out


def _quantile_edges_from_sorted(
    vals_sorted: list[float], bins: int, plow: float, phigh: float
) -> list[float]:
    if not vals_sorted:
        return []
    qs = [plow + (phigh - plow) * i / bins for i in range(bins + 1)]
    edges = [_quantile(vals_sorted, q) for q in qs]
    return _normalize_edges(edges)


def _category_quantile_edges(
    vals_sorted: list[float],
    bins: int,
    plow: float,
    phigh: float,
    global_lo: float,
    global_hi: float,
) -> list[float]:
    if not vals_sorted:
        return []
    edges = _quantile_edges_from_sorted(vals_sorted, bins, plow, phigh)
    if not edges:
        return []
    edges[0] = float(global_lo)
    edges[-1] = float(global_hi)
    return _normalize_edges(edges)


def _bin_idx(edges: list[float], val: float) -> int:
    if not edges or val < edges[0] or val > edges[-1]:
        return -1
    idx = bisect.bisect_right(edges, val) - 1
    if idx < 0:
        return -1
    if idx >= len(edges) - 1:
        return len(edges) - 2
    return idx


def _sample_from_hist(
    edges: list[float],
    counts: list[int],
    sample_size: int,
    rng: "np.random.Generator",
) -> list[float]:
    total = int(sum(counts))
    if total <= 0 or not edges or len(edges) < 2:
        return []
    size = sample_size if sample_size > 0 else total
    size = min(size, total)
    import numpy as np

    probs = np.asarray(counts, dtype=float)
    probs_sum = float(probs.sum())
    if probs_sum <= 0:
        return []
    probs = probs / probs_sum
    idxs = rng.choice(len(counts), size=size, replace=True, p=probs)
    left = np.asarray(edges[:-1], dtype=float)[idxs]
    right = np.asarray(edges[1:], dtype=float)[idxs]
    u = rng.random(size)
    return (left + (right - left) * u).tolist()


def plot_values_from_hist(
    hist_counts: dict[str, dict[str, list[int]]],
    hist_edges: dict[str, dict[str, list[float]]],
    scalar_name: str,
    sample_size: int,
    seed: int = 0,
) -> dict[str, list[float]]:
    try:
        import numpy as np
    except Exception:
        return {}
    rng = np.random.default_rng(seed)
    out: dict[str, list[float]] = {}
    for cat, name_map in hist_counts.items():
        counts = name_map.get(scalar_name)
        edges = hist_edges.get(cat, {}).get(scalar_name)
        if not counts or not edges:
            continue
        vals = _sample_from_hist(edges, counts, sample_size, rng)
        if vals:
            out[cat] = vals
    return out


def plot_values_from_hist_csvs(
    out_dir: Path,
    output_prefix: str,
    scalar_name: str,
    categories: Iterable[str],
    sample_size: int,
    seed: int = 0,
) -> dict[str, list[float]]:
    try:
        import numpy as np
    except Exception:
        return {}
    rng = np.random.default_rng(seed)
    out: dict[str, list[float]] = {}
    for cat in categories:
        csv_path = out_dir / f"{output_prefix}_{cat}_{scalar_name}_hist.csv"
        if not csv_path.exists():
            continue
        edges: list[float] = []
        counts: list[int] = []
        with open(csv_path, newline="", encoding="ascii") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                try:
                    left = float(row.get("bin_left", "nan"))
                    right = float(row.get("bin_right", "nan"))
                    cnt = int(float(row.get("count", "0")))
                except ValueError:
                    continue
                if not edges:
                    edges.append(left)
                edges.append(right)
                counts.append(cnt)
        if not edges or not counts:
            continue
        vals = _sample_from_hist(edges, counts, sample_size, rng)
        if vals:
            out[cat] = vals
    return out


def plot_values_from_raw_files(
    cache_dir: Path,
    output_prefix: str,
    scalar_name: str,
    categories: Iterable[str],
) -> dict[str, "np.ndarray"]:
    try:
        import numpy as np
    except Exception:
        return {}
    out: dict[str, "np.ndarray"] = {}
    for cat in categories:
        bin_path = cache_dir / f"{output_prefix}_{cat}_{scalar_name}.bin"
        if not bin_path.exists():
            continue
        size = bin_path.stat().st_size
        if size <= 0:
            continue
        count = size // 8
        if count <= 0:
            continue
        out[cat] = np.memmap(
            bin_path, dtype=np.float64, mode="r", shape=(count,)
        )
    return out


def read_inputs(
    paths: Iterable[Path], log10_field_value: bool = False
) -> tuple[list[str], dict[str, dict[str, list[float]]]]:
    scalars: list[str] = []
    hist_scalars: list[str] = []
    values: dict[str, dict[str, list[float]]] = {}
    for path in paths:
        with open(path, newline="", encoding="ascii") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames:
                continue
            if not scalars:
                scalar_candidates = [
                    name
                    for name in reader.fieldnames
                    if name
                    not in (
                        "delaunay_id",
                        "is_wall",
                        "is_filament",
                        "is_filament_manifold",
                        "is_cluster",
                        "is_cluster_manifold",
                    )
                ]
                scalars = scalar_candidates
                if "field_value" in scalars and "log10_field_value" not in scalars:
                    scalars.append("log10_field_value")
            cluster_key = None
            if "is_cluster" in reader.fieldnames:
                cluster_key = "is_cluster"
            elif "is_cluster_manifold" in reader.fieldnames:
                cluster_key = "is_cluster_manifold"
            for row in reader:
                is_wall = row.get("is_wall", "0").strip() == "1"
                is_filman = row.get("is_filament_manifold", "0").strip() == "1"
                is_cluster = (
                    cluster_key is not None and row.get(cluster_key, "0").strip() == "1"
                )

                include_filman = "is_filament_manifold" in reader.fieldnames
                include_cluster = cluster_key is not None
                categories = _categories(is_wall, is_filman, include_filman, is_cluster, include_cluster)

                field_val = None
                if "field_value" in row:
                    try:
                        field_val = float(row["field_value"])
                    except ValueError:
                        field_val = None
                for cat in categories:
                    values.setdefault(cat, {})
                    for name in scalars:
                        if name == "log10_field_value":
                            if field_val is None:
                                continue
                            val = _maybe_log10_field_value(field_val, True)
                            if val is None:
                                continue
                        else:
                            if name not in row:
                                continue
                            try:
                                val = float(row[name])
                            except ValueError:
                                continue
                        values[cat].setdefault(name, []).append(val)
    return scalars, values


def _collect(lazy_frame: "pl.LazyFrame") -> "pl.DataFrame":
    try:
        return lazy_frame.collect(engine="streaming")
    except TypeError:
        return lazy_frame.collect(streaming=True)


def polars_category_frame(lazy, include_filman: bool, include_cluster: bool) -> "pl.LazyFrame":
    import polars as pl

    is_wall = pl.col("is_wall") == 1
    is_filman = pl.col("is_filament_manifold") == 1 if include_filman else pl.lit(False)
    is_cluster = pl.col("is_cluster") == 1 if include_cluster else pl.lit(False)

    cat_exprs = [
        # Base categories
        pl.when(is_wall).then(pl.lit("walls")),
        pl.when(is_wall & ~is_filman & ~is_cluster).then(pl.lit("walls_only")),
        pl.when(is_filman).then(pl.lit("filmans")),
        pl.when(is_filman & ~is_wall & ~is_cluster).then(pl.lit("filmans_only")),
        pl.when(is_cluster).then(pl.lit("clusters")),
        pl.when(is_cluster & ~is_wall & ~is_filman).then(pl.lit("clusters_only")),
        pl.when(~is_wall & ~is_filman & ~is_cluster).then(pl.lit("unassigned")),
        # Exclusions
        pl.when(is_wall & ~is_cluster).then(pl.lit("walls_not_clusters")),
        pl.when(is_wall & ~is_filman).then(pl.lit("walls_not_filmans")),
        pl.when(is_filman & ~is_cluster).then(pl.lit("filmans_not_clusters")),
        pl.when(is_filman & ~is_wall).then(pl.lit("filmans_not_walls")),
        pl.when(is_cluster & ~is_filman).then(pl.lit("clusters_not_filmans")),
        pl.when(is_cluster & ~is_wall).then(pl.lit("clusters_not_walls")),
        # Universe
        pl.when(pl.lit(True)).then(pl.lit("unassigned_walls_filmans_clusters")),
        # Exclusive pairwise intersections (sigma-algebra atoms)
        pl.when(is_wall & is_filman & ~is_cluster).then(pl.lit("shared_walls_filmans")),
        pl.when(is_wall & is_cluster & ~is_filman).then(pl.lit("shared_walls_clusters")),
        pl.when(is_filman & is_cluster & ~is_wall).then(pl.lit("shared_filmans_clusters")),
        pl.when(is_wall & is_filman & is_cluster).then(pl.lit("shared_walls_filmans_clusters")),
        # Sigma-algebra: inclusive intersections
        pl.when(is_wall & is_filman).then(pl.lit("walls_and_filmans")),
        pl.when(is_wall & is_cluster).then(pl.lit("walls_and_clusters")),
        pl.when(is_filman & is_cluster).then(pl.lit("filmans_and_clusters")),
        # Pairwise unions
        pl.when(is_wall | is_filman).then(pl.lit("walls_or_filmans")),
        pl.when(is_wall | is_cluster).then(pl.lit("walls_or_clusters")),
        pl.when(is_filman | is_cluster).then(pl.lit("filmans_or_clusters")),
        pl.when(is_wall | is_filman | is_cluster).then(pl.lit("all_structures")),
        # Complements of generators
        pl.when(~is_wall).then(pl.lit("not_walls")),
        pl.when(~is_filman).then(pl.lit("not_filmans")),
        pl.when(~is_cluster).then(pl.lit("not_clusters")),
        # Complements of pairwise intersections
        pl.when(~(is_wall & is_filman)).then(pl.lit("not_walls_and_filmans")),
        pl.when(~(is_wall & is_cluster)).then(pl.lit("not_walls_and_clusters")),
        pl.when(~(is_filman & is_cluster)).then(pl.lit("not_filmans_and_clusters")),
        # Complement of triple intersection
        pl.when(~(is_wall & is_filman & is_cluster)).then(pl.lit("not_walls_and_filmans_and_clusters")),
        # Complements of pairwise unions
        pl.when(~is_wall & ~is_filman).then(pl.lit("neither_walls_nor_filmans")),
        pl.when(~is_wall & ~is_cluster).then(pl.lit("neither_walls_nor_clusters")),
        pl.when(~is_filman & ~is_cluster).then(pl.lit("neither_filmans_nor_clusters")),
    ]
    category_col = pl.concat_list(cat_exprs).list.drop_nulls().alias("category")
    return lazy.with_columns(category_col).explode("category")


def _prepare_polars_lazy(lazy: "pl.LazyFrame", include_filman: bool, include_cluster: bool) -> "pl.LazyFrame":
    import polars as pl
    cols = set(lazy.collect_schema().names())
    if "field_value" in cols and "log10_field_value" not in cols:
        lazy = lazy.with_columns(
            pl.when(pl.col("field_value") > 0)
            .then(pl.col("field_value").log10())
            .otherwise(None)
            .alias("log10_field_value")
        )
        cols = set(lazy.collect_schema().names())
    if not include_filman:
        lazy = lazy.with_columns(pl.lit(0).alias("is_filament_manifold"))
    if "is_cluster" not in cols:
        if "is_cluster_manifold" in cols:
            lazy = lazy.with_columns(pl.col("is_cluster_manifold").alias("is_cluster"))
        else:
            lazy = lazy.with_columns(pl.lit(0).alias("is_cluster"))
    return lazy


def _check_matplotlib(no_plots: bool) -> bool:
    if no_plots:
        return False
    try:
        import matplotlib.pyplot as plt  # noqa: F401
        return True
    except Exception:
        return False


def _write_composition_and_transform_stats(
    out_dir: Path,
    prefix: str,
    stats_rows: list[dict],
    scalars: list[str],
    have_plots: bool,
    args: argparse.Namespace,
) -> None:
    _write_transform_stats(out_dir, prefix, stats_rows, scalars)
    composition_counts = _composition_counts_from_rows(stats_rows)
    composition_sums = _composition_sums_from_rows(stats_rows, "field_value")
    write_filman_walls_composition_plot(
        out_dir, prefix, composition_counts, have_plots, args.plot_fontscale, args.plot_dpi,
    )
    write_filman_walls_composition_mass_plot(
        out_dir, prefix, composition_sums, have_plots, args.plot_fontscale, args.plot_dpi,
        scalar_label="field_value",
    )


def _write_hist_plots(
    out_dir: Path,
    prefix: str,
    hist_counts: dict[str, dict[str, list[int]]],
    hist_edges: dict[str, dict[str, list[float]]],
    medians_map: dict[str, dict[str, float]],
    cat_counts: dict[str, int],
    have_plots: bool,
    args: argparse.Namespace,
) -> None:
    for cat, name_map in hist_counts.items():
        for name, hist in name_map.items():
            edges = hist_edges.get(cat, {}).get(name)
            if not edges:
                continue
            csv_path = out_dir / f"{prefix}_{cat}_{name}_hist.csv"
            write_histogram_csv(csv_path, edges, hist, cat_counts.get(cat, 0))
            if have_plots:
                import matplotlib.pyplot as plt
                centers = [(edges[i] + edges[i + 1]) / 2 for i in range(len(edges) - 1)]
                widths = [edges[i + 1] - edges[i] for i in range(len(edges) - 1)]
                fig = plt.figure()
                plt.bar(centers, hist, width=widths, align="center")
                median_val = medians_map.get(cat, {}).get(name)
                if median_val is not None:
                    plt.axvline(median_val, color="red", linestyle="--", linewidth=1)
                plt.xlabel(_scalar_label(name))
                plt.ylabel("count")
                plt.title(f"{cat}: {_scalar_label(name)} (n={cat_counts.get(cat, 0):.2e})")
                fig_path = out_dir / f"{prefix}_{cat}_{name}_hist.png"
                fig.savefig(fig_path, dpi=args.plot_dpi, bbox_inches="tight")
                plt.close(fig)


def run_polars(paths: list[Path], args: argparse.Namespace) -> None:
    import polars as pl
    violin_scalar = _plot_scalar_name(args.violin_scalar, args.log10_field_value)
    box_scalar = _plot_scalar_name(args.box_scalar, args.log10_field_value)

    scalars, include_filman, include_cluster = _infer_schema_from_headers(paths)
    if not scalars:
        raise SystemExit("[error] No scalar columns found in topology_points.csv inputs.")
    lazy = _prepare_polars_lazy(pl.concat([pl.scan_csv(str(p)) for p in paths]), include_filman, include_cluster)

    data = polars_category_frame(lazy, include_filman, include_cluster)
    agg_exprs = [pl.len().alias("count")]
    for name in scalars:
        agg_exprs.extend(
            [
                pl.col(name).sum().alias(f"{name}_sum"),
                pl.col(name).mean().alias(f"{name}_mean"),
                pl.col(name).min().alias(f"{name}_min"),
                pl.col(name).quantile(0.25, "linear").alias(f"{name}_q25"),
                pl.col(name).quantile(0.5, "linear").alias(f"{name}_median"),
                pl.col(name).quantile(0.75, "linear").alias(f"{name}_q75"),
                pl.col(name).max().alias(f"{name}_max"),
                pl.col(name).std(ddof=0).alias(f"{name}_std"),
            ]
        )
    stats_df = _collect(data.group_by("category").agg(agg_exprs))
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    stats_path = out_dir / f"{args.output_prefix}_topology_stats.csv"
    fieldnames = _build_fieldnames(scalars)
    stats_df = stats_df.select(fieldnames)
    stats_rows = stats_df.to_dicts()
    cat_counts = {row["category"]: int(row["count"]) for row in stats_rows}
    stats_rows = _seed_category_rows(stats_rows, scalars, ALL_CATEGORIES)
    stats_df = pl.DataFrame(stats_rows, schema=fieldnames)
    stats_df.write_csv(stats_path)
    print(f"[done] wrote {stats_path}")
    have_plots = _check_matplotlib(args.no_plots)
    _write_composition_and_transform_stats(out_dir, args.output_prefix, stats_rows, scalars, have_plots, args)
    medians_map: dict[str, dict[str, float]] = {}
    box_stats_map: dict[str, dict[str, dict[str, float]]] = {}
    for row in stats_rows:
        cat = row["category"]
        medians_map[cat] = {}
        for name in scalars:
            medians_map[cat][name] = float(row.get(f"{name}_median", 0.0) or 0.0)
            box_stats_map.setdefault(name, {})[cat] = {
                "count": float(row.get("count", 0.0) or 0.0),
                "mean": float(row.get(f"{name}_mean", 0.0) or 0.0),
                "min": float(row.get(f"{name}_min", 0.0) or 0.0),
                "q25": float(row.get(f"{name}_q25", 0.0) or 0.0),
                "median": float(row.get(f"{name}_median", 0.0) or 0.0),
                "q75": float(row.get(f"{name}_q75", 0.0) or 0.0),
                "max": float(row.get(f"{name}_max", 0.0) or 0.0),
                "std": float(row.get(f"{name}_std", 0.0) or 0.0),
            }

    hist_scalars = _select_hist_scalars(scalars, args)
    plow = phigh = None
    if args.hist_percentile_range:
        plow, phigh = args.hist_percentile_range
        plow = plow / 100.0
        phigh = phigh / 100.0
    qlo = 0.0 if plow is None else plow
    qhi = 1.0 if phigh is None else phigh

    hist_counts: dict[str, dict[str, list[int]]] = {}
    hist_edges: dict[str, dict[str, list[float]]] = {}
    plot_values_all: dict[str, dict[str, list[float]]] = {}

    if args.hist_bin_mode == "global":
        try:
            import numpy as np
        except Exception as exc:
            raise SystemExit(f"[error] numpy is required for histograms: {exc}")

        for name in hist_scalars:
            if plow is not None and phigh is not None:
                row = _collect(
                    lazy.select(
                        [
                            pl.col(name).quantile(plow, "linear").alias("min"),
                            pl.col(name).quantile(phigh, "linear").alias("max"),
                        ]
                    )
                ).to_dicts()[0]
            else:
                row = _collect(
                    lazy.select(
                        [pl.col(name).min().alias("min"), pl.col(name).max().alias("max")]
                    )
                ).to_dicts()[0]
            if row.get("min") is None or row.get("max") is None:
                continue
            low = float(row["min"])
            high = float(row["max"])
            if high <= low:
                continue
            edges = np.linspace(low, high, args.bins + 1, dtype=float).tolist()
            edges_arr = np.asarray(edges, dtype=float)

            def _bin_series(series: "pl.Series") -> "pl.Series":
                arr = series.to_numpy()
                idx = np.searchsorted(edges_arr, arr, side="right") - 1
                idx[(arr < edges_arr[0]) | (arr > edges_arr[-1])] = -1
                return pl.Series(idx, dtype=pl.Int64)

            bin_expr = pl.col(name).map_batches(_bin_series, return_dtype=pl.Int64).alias("bin")
            grouped = _collect(
                data.select(["category", bin_expr])
                .filter(pl.col("bin") >= 0)
                .group_by(["category", "bin"])
                .agg(pl.len().alias("count"))
            )

            for cat in cat_counts:
                counts = [0] * (len(edges) - 1)
                hist_counts.setdefault(cat, {})[name] = counts
                hist_edges.setdefault(cat, {})[name] = edges
            for row in grouped.to_dicts():
                cat = row["category"]
                bin_idx = int(row["bin"])
                cnt = int(row["count"])
                if cat in hist_counts and name in hist_counts[cat]:
                    if 0 <= bin_idx < len(hist_counts[cat][name]):
                        hist_counts[cat][name][bin_idx] = cnt

        _write_hist_plots(out_dir, args.output_prefix, hist_counts, hist_edges, medians_map, cat_counts, have_plots, args)
    else:
        # Per-category quantile bins with global support per scalar.
        global_support: dict[str, tuple[float, float]] = {}
        for name in hist_scalars:
            if plow is not None and phigh is not None:
                row = _collect(
                    lazy.select(
                        [
                            pl.col(name).quantile(plow, "linear").alias("min"),
                            pl.col(name).quantile(phigh, "linear").alias("max"),
                        ]
                    )
                ).to_dicts()[0]
            else:
                row = _collect(
                    lazy.select(
                        [pl.col(name).min().alias("min"), pl.col(name).max().alias("max")]
                    )
                ).to_dicts()[0]
            if row.get("min") is None or row.get("max") is None:
                continue
            global_support[name] = (float(row["min"]), float(row["max"]))

        plot_values_all = collect_plot_values(paths, scalars, args.log10_field_value)
        for name in hist_scalars:
            if name not in global_support:
                continue
            gmin, gmax = global_support[name]
            for cat in cat_counts:
                vals = plot_values_all.get(name, {}).get(cat, [])
                if not vals:
                    continue
                vals_sorted = sorted(vals)
                edges = _category_quantile_edges(vals_sorted, args.bins, qlo, qhi, gmin, gmax)
                if not edges:
                    continue
                counts = [0] * (len(edges) - 1)
                for v in vals:
                    idx = _bin_idx(edges, v)
                    if idx >= 0:
                        counts[idx] += 1
                hist_counts.setdefault(cat, {})[name] = counts
                hist_edges.setdefault(cat, {})[name] = edges
        _write_hist_plots(out_dir, args.output_prefix, hist_counts, hist_edges, medians_map, cat_counts, have_plots, args)

    if have_plots:
        plot_scalars = {violin_scalar, box_scalar}
        plot_scalars = {s for s in plot_scalars if s in scalars}
        if plot_scalars:
            plot_values: dict[str, dict[str, list[float]]] = {}
            if not args.plots_from_raw and hist_counts and hist_edges:
                for name in plot_scalars:
                    if name in hist_scalars:
                        plot_values[name] = plot_values_from_hist(
                            hist_counts, hist_edges, name, args.plot_sample_size
                        )
                    else:
                        if not plot_values_all:
                            plot_values_all = collect_plot_values(paths, list(plot_scalars), args.log10_field_value)
                        plot_values[name] = plot_values_all.get(name, {})
            else:
                if not plot_values_all:
                    plot_values_all = collect_plot_values(paths, list(plot_scalars), args.log10_field_value)
                plot_values = {name: plot_values_all.get(name, {}) for name in plot_scalars}
            if violin_scalar in plot_values:
                write_filman_walls_violin_plot(
                    out_dir,
                    args.output_prefix,
                    plot_values[violin_scalar],
                    have_plots,
                    scalar_name=violin_scalar,
                    font_scale=args.plot_fontscale,
                    plot_dpi=args.plot_dpi,
                    plot_percentile_range=tuple(args.plot_percentile_range)
                    if args.plot_percentile_range
                    else None,
                    log10_field_value=args.log10_field_value,
                )
            if box_scalar in scalars:
                if args.plot_percentile_range and box_scalar in plot_values:
                    exact_box_stats: dict[str, dict[str, float]] = {}
                    keys, _, _ = _composition_lists()
                    for key in keys:
                        vals = plot_values[box_scalar].get(key, [])
                        vals = trim_by_percentile(
                            vals, args.plot_percentile_range[0], args.plot_percentile_range[1]
                        )
                        exact_box_stats[key] = aggregate_values_np(vals)
                    write_filman_walls_box_plot(
                        out_dir,
                        args.output_prefix,
                        exact_box_stats,
                        have_plots,
                        scalar_name=box_scalar,
                        font_scale=args.plot_fontscale,
                        plot_dpi=args.plot_dpi,
                        log10_field_value=args.log10_field_value,
                    )
                else:
                    write_filman_walls_box_plot(
                        out_dir,
                        args.output_prefix,
                        box_stats_map.get(box_scalar, {}),
                        have_plots,
                        scalar_name=box_scalar,
                        font_scale=args.plot_fontscale,
                        plot_dpi=args.plot_dpi,
                        log10_field_value=args.log10_field_value,
                    )


def _categories(
    is_wall: bool,
    is_filman: bool,
    include_filman: bool,
    is_cluster: bool,
    include_cluster: bool,
) -> list[str]:
    cats: list[str] = []
    # Base categories
    if is_wall:
        cats.append("walls")
    if is_wall and not is_filman and not is_cluster:
        cats.append("walls_only")
    if is_filman:
        cats.append("filmans")
    if is_filman and not is_wall and not is_cluster:
        cats.append("filmans_only")
    if is_cluster:
        cats.append("clusters")
    if is_cluster and not is_wall and not is_filman:
        cats.append("clusters_only")
    if not is_wall and not is_filman and not is_cluster:
        cats.append("unassigned")
    # Exclusions
    if is_wall and not is_cluster:
        cats.append("walls_not_clusters")
    if is_wall and not is_filman:
        cats.append("walls_not_filmans")
    if is_filman and not is_cluster:
        cats.append("filmans_not_clusters")
    if is_filman and not is_wall:
        cats.append("filmans_not_walls")
    if is_cluster and not is_filman:
        cats.append("clusters_not_filmans")
    if is_cluster and not is_wall:
        cats.append("clusters_not_walls")
    # Universe (all points — unconditional after removing arc filament dimension)
    cats.append("unassigned_walls_filmans_clusters")
    # Exclusive pairwise intersections (atoms of the sigma-algebra)
    if is_wall and is_filman and not is_cluster:
        cats.append("shared_walls_filmans")
    if is_wall and is_cluster and not is_filman:
        cats.append("shared_walls_clusters")
    if is_filman and is_cluster and not is_wall:
        cats.append("shared_filmans_clusters")
    if is_wall and is_filman and is_cluster:
        cats.append("shared_walls_filmans_clusters")
    # Sigma-algebra generated by {W, F, C}
    # Inclusive pairwise intersections
    if is_wall and is_filman:
        cats.append("walls_and_filmans")
    if is_wall and is_cluster:
        cats.append("walls_and_clusters")
    if is_filman and is_cluster:
        cats.append("filmans_and_clusters")
    # Pairwise unions
    if is_wall or is_filman:
        cats.append("walls_or_filmans")
    if is_wall or is_cluster:
        cats.append("walls_or_clusters")
    if is_filman or is_cluster:
        cats.append("filmans_or_clusters")
    if is_wall or is_filman or is_cluster:
        cats.append("all_structures")
    # Complements of generators
    if not is_wall:
        cats.append("not_walls")
    if not is_filman:
        cats.append("not_filmans")
    if not is_cluster:
        cats.append("not_clusters")
    # Complements of pairwise intersections
    if not (is_wall and is_filman):
        cats.append("not_walls_and_filmans")
    if not (is_wall and is_cluster):
        cats.append("not_walls_and_clusters")
    if not (is_filman and is_cluster):
        cats.append("not_filmans_and_clusters")
    # Complement of triple intersection
    if not (is_wall and is_filman and is_cluster):
        cats.append("not_walls_and_filmans_and_clusters")
    # Complements of pairwise unions
    if not is_wall and not is_filman:
        cats.append("neither_walls_nor_filmans")
    if not is_wall and not is_cluster:
        cats.append("neither_walls_nor_clusters")
    if not is_filman and not is_cluster:
        cats.append("neither_filmans_nor_clusters")
    return cats


def collect_plot_values(
    paths: Iterable[Path], scalars_needed: list[str], log10_field_value: bool = False
) -> dict[str, dict[str, list[float]]]:
    values: dict[str, dict[str, list[float]]] = {name: {} for name in scalars_needed}
    scalars_set = set(scalars_needed)
    for path in paths:
        with open(path, newline="", encoding="ascii") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames:
                continue
            include_filman = "is_filament_manifold" in reader.fieldnames
            include_cluster = "is_cluster" in reader.fieldnames or "is_cluster_manifold" in reader.fieldnames
            for row in reader:
                is_wall = row.get("is_wall", "0").strip() == "1"
                is_filman = (
                    include_filman and row.get("is_filament_manifold", "0").strip() == "1"
                )
                cluster_key = "is_cluster" if "is_cluster" in reader.fieldnames else "is_cluster_manifold"
                is_cluster = include_cluster and row.get(cluster_key, "0").strip() == "1"
                cats = _categories(is_wall, is_filman, include_filman, is_cluster, include_cluster)
                field_val = None
                if "field_value" in row:
                    try:
                        field_val = float(row["field_value"])
                    except ValueError:
                        field_val = None
                for name in scalars_set:
                    if name == "log10_field_value":
                        if field_val is None:
                            continue
                        val = _maybe_log10_field_value(field_val, True)
                        if val is None:
                            continue
                    else:
                        if name not in row:
                            continue
                        try:
                            val = float(row[name])
                        except ValueError:
                            continue
                    cat_map = values.setdefault(name, {})
                    for cat in cats:
                        cat_map.setdefault(cat, []).append(val)
    return values


def run_stream(paths: list[Path], args: argparse.Namespace) -> None:
    import random
    violin_scalar = _plot_scalar_name(args.violin_scalar, args.log10_field_value)
    box_scalar = _plot_scalar_name(args.box_scalar, args.log10_field_value)

    sample_size = 20000
    scalars, include_filman, include_cluster = _infer_schema_from_headers(paths)
    if not scalars:
        raise SystemExit("[error] No scalar columns found in topology_points.csv inputs.")
    hist_scalars = _select_hist_scalars(scalars, args)
    stats: dict[str, dict[str, dict[str, float]]] = {}
    category_count: dict[str, int] = {}
    samples: dict[str, dict[str, list[float]]] = {}
    global_samples: dict[str, list[float]] = {}
    global_counts: dict[str, int] = {}
    def _stat_store(cat: str, name: str) -> dict[str, float]:
        return stats.setdefault(cat, {}).setdefault(
            name,
            {
                "count": 0.0,
                "sum": 0.0,
                "mean": 0.0,
                "m2": 0.0,
                "min": float("inf"),
                "max": float("-inf"),
            },
        )

    def _sample_store(cat: str, name: str) -> list[float]:
        return samples.setdefault(cat, {}).setdefault(name, [])

    def _global_sample_store(name: str) -> list[float]:
        return global_samples.setdefault(name, [])

    # pass 1: streaming stats + samples
    for path in paths:
        with open(path, newline="", encoding="ascii") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames:
                continue
            cluster_key = "is_cluster" if "is_cluster" in reader.fieldnames else "is_cluster_manifold"
            for row in reader:
                is_wall = row.get("is_wall", "0").strip() == "1"
                is_filman = include_filman and row.get("is_filament_manifold", "0").strip() == "1"
                is_cluster = include_cluster and row.get(cluster_key, "0").strip() == "1"
                cats = _categories(is_wall, is_filman, include_filman, is_cluster, include_cluster)
                field_val = None
                if "field_value" in row:
                    try:
                        field_val = float(row.get("field_value", "0") or 0.0)
                    except ValueError:
                        field_val = None
                for cat in cats:
                    category_count[cat] = category_count.get(cat, 0) + 1
                for name in scalars:
                    if name == "log10_field_value":
                        if field_val is None:
                            continue
                        val = _maybe_log10_field_value(field_val, True)
                        if val is None:
                            continue
                    else:
                        try:
                            val = float(row.get(name, "0") or 0.0)
                        except ValueError:
                            continue
                    for cat in cats:
                        store = _stat_store(cat, name)
                        store["count"] += 1.0
                        store["sum"] += val
                        if val < store["min"]:
                            store["min"] = val
                        if val > store["max"]:
                            store["max"] = val
                        delta = val - store["mean"]
                        store["mean"] += delta / store["count"]
                        delta2 = val - store["mean"]
                        store["m2"] += delta * delta2
                        samp = _sample_store(cat, name)
                        seen = int(store["count"])
                        if len(samp) < sample_size:
                            samp.append(val)
                        else:
                            j = random.randint(1, seen)
                            if j <= sample_size:
                                samp[j - 1] = val
                        # global reservoir sample for quantile bins
                        if name in hist_scalars:
                            gs = _global_sample_store(name)
                            global_counts[name] = global_counts.get(name, 0) + 1
                            seen_global = global_counts[name]
                            if len(gs) < sample_size:
                                gs.append(val)
                            else:
                                j = random.randint(1, seen_global)
                                if j <= sample_size:
                                    gs[j - 1] = val

    if not scalars:
        raise SystemExit("[error] No scalar columns found in topology_points.csv inputs.")

    # quantile bin edges per category, with global support per scalar
    plow = phigh = None
    if args.hist_percentile_range:
        plow, phigh = args.hist_percentile_range
        plow = plow / 100.0
        phigh = phigh / 100.0
    qlo = 0.0 if plow is None else plow
    qhi = 1.0 if phigh is None else phigh

    global_support: dict[str, tuple[float, float]] = {}
    for name in hist_scalars:
        low = float("inf")
        high = float("-inf")
        for scal_map in stats.values():
            store = scal_map.get(name)
            if not store or store["count"] <= 0:
                continue
            low = min(low, store["min"])
            high = max(high, store["max"])
        if low == float("inf") or high == float("-inf") or high <= low:
            continue
        if plow is not None and phigh is not None:
            gs = sorted(global_samples.get(name, []))
            if gs:
                low = _quantile(gs, plow)
                high = _quantile(gs, phigh)
        if high > low:
            global_support[name] = (low, high)


    edges_by_cat: dict[str, dict[str, list[float]]] = {}
    if args.hist_bin_mode == "global":
        for name in hist_scalars:
            if name not in global_support:
                continue
            gmin, gmax = global_support[name]
            if gmax <= gmin:
                continue
            edges = [gmin + (gmax - gmin) * i / args.bins for i in range(args.bins + 1)]
            for cat in category_count:
                edges_by_cat.setdefault(cat, {})[name] = edges
    else:
        for cat in category_count:
            for name in hist_scalars:
                if name not in global_support:
                    continue
                samples_cat = sorted(samples.get(cat, {}).get(name, []))
                if not samples_cat:
                    continue
                gmin, gmax = global_support[name]
                edges = _category_quantile_edges(samples_cat, args.bins, qlo, qhi, gmin, gmax)
                if edges:
                    edges_by_cat.setdefault(cat, {})[name] = edges

    # pass 2: histograms
    counts: dict[str, dict[str, list[int]]] = {}
    for path in paths:
        with open(path, newline="", encoding="ascii") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames:
                continue
            cluster_key = "is_cluster" if "is_cluster" in reader.fieldnames else "is_cluster_manifold"
            for row in reader:
                is_wall = row.get("is_wall", "0").strip() == "1"
                is_filman = include_filman and row.get("is_filament_manifold", "0").strip() == "1"
                is_cluster = include_cluster and row.get(cluster_key, "0").strip() == "1"
                cats = _categories(is_wall, is_filman, include_filman, is_cluster, include_cluster)
                field_val = None
                if "field_value" in row:
                    try:
                        field_val = float(row.get("field_value", "0") or 0.0)
                    except ValueError:
                        field_val = None
                for name in scalars:
                    try:
                        if name == "log10_field_value":
                            if field_val is None:
                                continue
                            val = _maybe_log10_field_value(field_val, True)
                            if val is None:
                                continue
                        else:
                            val = float(row.get(name, "0") or 0.0)
                    except ValueError:
                        continue
                    if name not in hist_scalars:
                        continue
                    for cat in cats:
                        edges = edges_by_cat.get(cat, {}).get(name)
                        if not edges:
                            continue
                        if val < edges[0] or val > edges[-1]:
                            continue
                        bin_idx = bisect.bisect_right(edges, val) - 1
                        if bin_idx < 0 or bin_idx >= len(edges) - 1:
                            continue
                        counts.setdefault(cat, {}).setdefault(name, [0] * (len(edges) - 1))
                        counts[cat][name][bin_idx] += 1


    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    stats_path = out_dir / f"{args.output_prefix}_topology_stats.csv"
    fieldnames = _build_fieldnames(scalars)
    stats_rows: list[dict[str, object]] = []
    with open(stats_path, "w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        medians_map: dict[str, dict[str, float]] = {}
        box_stats_map: dict[str, dict[str, dict[str, float]]] = {}
        for cat in ALL_CATEGORIES:
            row: dict[str, object] = {"category": cat, "count": category_count.get(cat, 0)}
            for name in scalars:
                store = stats.get(cat, {}).get(name)
                if not store or store["count"] <= 0:
                    for key in ("sum", "mean", "min", "q25", "median", "q75", "max", "std"):
                        row[f"{name}_{key}"] = 0.0
                    continue
                row[f"{name}_sum"] = store["sum"]
                row[f"{name}_mean"] = store["mean"]
                row[f"{name}_min"] = store["min"]
                row[f"{name}_max"] = store["max"]
                row[f"{name}_std"] = math.sqrt(store["m2"] / store["count"]) if store["count"] else 0.0
                qvals = {"q25": 0.0, "median": 0.0, "q75": 0.0}
                samp = samples.get(cat, {}).get(name, [])
                if samp:
                    samp_sorted = sorted(samp)
                    qvals["q25"] = _quantile(samp_sorted, 0.25)
                    qvals["median"] = _quantile(samp_sorted, 0.5)
                    qvals["q75"] = _quantile(samp_sorted, 0.75)
                row[f"{name}_q25"] = qvals["q25"]
                row[f"{name}_median"] = qvals["median"]
                row[f"{name}_q75"] = qvals["q75"]
                medians_map.setdefault(cat, {})[name] = qvals["median"]
                box_stats_map.setdefault(name, {})[cat] = {
                    "count": float(store["count"]),
                    "mean": float(store["mean"]),
                    "min": float(store["min"]),
                    "q25": float(qvals["q25"]),
                    "median": float(qvals["median"]),
                    "q75": float(qvals["q75"]),
                    "max": float(store["max"]),
                    "std": float(math.sqrt(store["m2"] / store["count"])) if store["count"] else 0.0,
                }
            writer.writerow(row)
            stats_rows.append(row)
    print(f"[done] wrote {stats_path}")
    have_plots = _check_matplotlib(args.no_plots)
    cat_counts = category_count
    _write_composition_and_transform_stats(out_dir, args.output_prefix, stats_rows, scalars, have_plots, args)
    _write_hist_plots(out_dir, args.output_prefix, counts, edges_by_cat, medians_map, cat_counts, have_plots, args)
    if have_plots:
        plot_scalars = {violin_scalar, box_scalar}
        plot_scalars = {s for s in plot_scalars if s in scalars}
        if plot_scalars:
            if not args.plots_from_raw and counts and edges_by_cat:
                plot_values = {}
                for name in plot_scalars:
                    if name in hist_scalars:
                        plot_values[name] = plot_values_from_hist(
                            counts, edges_by_cat, name, args.plot_sample_size
                        )
                    else:
                        plot_values[name] = collect_plot_values(
                            paths, [name], args.log10_field_value
                        ).get(name, {})
            else:
                plot_values = collect_plot_values(
                    paths, sorted(plot_scalars), args.log10_field_value
                )
            if violin_scalar in plot_values:
                write_filman_walls_violin_plot(
                    out_dir,
                    args.output_prefix,
                    plot_values[violin_scalar],
                    have_plots,
                    scalar_name=violin_scalar,
                    font_scale=args.plot_fontscale,
                    plot_dpi=args.plot_dpi,
                    plot_percentile_range=tuple(args.plot_percentile_range)
                    if args.plot_percentile_range
                    else None,
                    log10_field_value=args.log10_field_value,
                )
            if box_scalar in scalars:
                if args.plot_percentile_range and box_scalar in plot_values:
                    exact_box_stats: dict[str, dict[str, float]] = {}
                    keys, _, _ = _composition_lists()
                    for key in keys:
                        vals = plot_values[box_scalar].get(key, [])
                        vals = trim_by_percentile(
                            vals, args.plot_percentile_range[0], args.plot_percentile_range[1]
                        )
                        exact_box_stats[key] = aggregate_values_np(vals)
                    write_filman_walls_box_plot(
                        out_dir,
                        args.output_prefix,
                        exact_box_stats,
                        have_plots,
                        scalar_name=box_scalar,
                        font_scale=args.plot_fontscale,
                        plot_dpi=args.plot_dpi,
                        log10_field_value=args.log10_field_value,
                    )
                else:
                    write_filman_walls_box_plot(
                        out_dir,
                        args.output_prefix,
                        box_stats_map.get(box_scalar, {}),
                        have_plots,
                        scalar_name=box_scalar,
                        font_scale=args.plot_fontscale,
                        plot_dpi=args.plot_dpi,
                        log10_field_value=args.log10_field_value,
                    )


def write_histogram_csv(
    path: Path,
    edges: list[float],
    counts: list[int],
    total_count: int,
) -> None:
    with open(path, "w", newline="", encoding="ascii") as handle:
        writer = csv.writer(handle)
        writer.writerow(["bin_left", "bin_right", "count", "total_count"])
        for i, count in enumerate(counts):
            writer.writerow([edges[i], edges[i + 1], count, total_count])


def _hist_quantile_from_counts(
    counts_list: list[int], edges_list: list[float], q: float
) -> float:
    total = int(sum(counts_list))
    if total <= 0:
        return float("nan")
    try:
        import numpy as np
    except Exception:
        return float("nan")
    counts_arr = np.array(counts_list, dtype=float)
    cdf = np.cumsum(counts_arr)
    target = q * total
    idx = int(np.searchsorted(cdf, target, side="left"))
    if idx < 0:
        idx = 0
    if idx >= len(counts_arr):
        idx = len(counts_arr) - 1
    count_in_bin = counts_arr[idx]
    if count_in_bin <= 0:
        return float((edges_list[idx] + edges_list[idx + 1]) / 2)
    prev = cdf[idx - 1] if idx > 0 else 0.0
    frac = (target - prev) / count_in_bin
    left = edges_list[idx]
    right = edges_list[idx + 1]
    return float(left + frac * (right - left))


def _hist_stats_from_counts(counts_list: list[int], edges_list: list[float]) -> dict[str, float]:
    total = int(sum(counts_list))
    if total <= 0:
        return {
            "count": 0,
            "mean": float("nan"),
            "min": float("nan"),
            "q25": float("nan"),
            "q50": float("nan"),
            "q75": float("nan"),
            "max": float("nan"),
        }
    try:
        import numpy as np
    except Exception:
        return {
            "count": total,
            "mean": float("nan"),
            "min": float(edges_list[0]),
            "q25": float("nan"),
            "q50": float("nan"),
            "q75": float("nan"),
            "max": float(edges_list[-1]),
        }

    centers = np.array(
        [(edges_list[i] + edges_list[i + 1]) / 2 for i in range(len(counts_list))],
        dtype=float,
    )
    counts_arr = np.array(counts_list, dtype=float)
    mean = float((counts_arr * centers).sum() / total)
    first_idx = next((i for i, c in enumerate(counts_list) if c > 0), 0)
    last_idx = len(counts_list) - 1 - next(
        (i for i, c in enumerate(reversed(counts_list)) if c > 0),
        0,
    )
    min_val = float(edges_list[first_idx])
    max_val = float(edges_list[last_idx + 1])
    return {
        "count": total,
        "mean": mean,
        "min": min_val,
        "q25": _hist_quantile_from_counts(counts_list, edges_list, 0.25),
        "q50": _hist_quantile_from_counts(counts_list, edges_list, 0.50),
        "q75": _hist_quantile_from_counts(counts_list, edges_list, 0.75),
        "max": max_val,
    }


def write_filman_walls_composition_plot(
    out_dir: Path,
    output_prefix: str,
    cat_counts: dict[str, int],
    have_plots: bool,
    font_scale: float,
    plot_dpi: int,
) -> None:
    if not have_plots:
        return
    keys, labels, colors = _composition_lists()
    total = sum(cat_counts.get(key, 0) for key in keys)
    if total <= 0:
        return
    percents = [cat_counts.get(key, 0) / total * 100.0 for key in keys]
    percents[-1] = 100.0 - sum(percents[:-1])

    import matplotlib.pyplot as plt

    sizes = _font_sizes(font_scale)
    fig, ax = plt.subplots(figsize=(9, 1.6))
    left = 0.0
    counts = []
    for key, label, color, pct in zip(keys, labels, colors, percents):
        bar = ax.barh([0], [pct], left=left, color=color, edgecolor="none", alpha=1.0)
        counts.append(cat_counts.get(key, 0))
        left += pct
    ax.set_xlim(0, 100)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_xlabel("")
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title("Particle membership in cosmic structures", fontsize=sizes["title"])
    table_bbox = [0.13, -0.8, 0.85, 0.6]
    table = ax.table(
        cellText=[
            [f"{c:.2e}" for c in counts],
            [f"{p:.2f}%" for p in percents],
        ],
        rowLabels=["count", "proportion"],
        colLabels=labels,
        cellLoc="center",
        rowLoc="center",
        loc="bottom",
        bbox=table_bbox,
    )
    table.auto_set_font_size(False)
    table.set_fontsize(sizes["table"])
    table.scale(0.9, 1.45)
    for idx, color in enumerate(colors):
        for row_idx in (0, 1, 2):
            cell = table[row_idx, idx]
            cell.get_text().set_color(color)
            _apply_text_outline(cell.get_text(), linewidth=0.3)
    for cell in table.get_celld().values():
        cell.set_linewidth(0.0)
        cell.set_edgecolor("none")
        cell.set_facecolor("none")
    ax.set_position([0.1, 0.52, 0.8, 0.4])
    fig.tight_layout()
    out_path = out_dir / f"{output_prefix}_filman_walls_composition.png"
    _set_transparent(fig, ax)
    fig.savefig(out_path, dpi=plot_dpi, bbox_inches="tight", transparent=True)
    plt.close(fig)


def write_filman_walls_composition_mass_plot(
    out_dir: Path,
    output_prefix: str,
    cat_sums: dict[str, float],
    have_plots: bool,
    font_scale: float,
    plot_dpi: int,
    scalar_label: str = "field_value",
) -> None:
    if not have_plots:
        return
    keys, labels, colors = _composition_lists()
    total = sum(cat_sums.get(key, 0.0) for key in keys)
    if total <= 0:
        return
    percents = [cat_sums.get(key, 0.0) / total * 100.0 for key in keys]
    percents[-1] = 100.0 - sum(percents[:-1])

    import matplotlib.pyplot as plt

    sizes = _font_sizes(font_scale)
    fig, ax = plt.subplots(figsize=(9, 1.6))
    left = 0.0
    sums = []
    for key, label, color, pct in zip(keys, labels, colors, percents):
        ax.barh([0], [pct], left=left, color=color, edgecolor="none", alpha=1.0)
        sums.append(cat_sums.get(key, 0.0))
        left += pct
    ax.set_xlim(0, 100)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_xlabel("")
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(
        f"Mass distribution across cosmic structures",
        fontsize=sizes["title"],
    )
    table_bbox = [0.13, -0.8, 0.85, 0.6]
    table = ax.table(
        cellText=[
            [f"{c:.2e}" for c in sums],
            [f"{p:.2f}%" for p in percents],
        ],
        rowLabels=["mass", "proportion"],
        colLabels=labels,
        cellLoc="center",
        rowLoc="center",
        loc="bottom",
        bbox=table_bbox,
    )
    table.auto_set_font_size(False)
    table.set_fontsize(sizes["table"])
    table.scale(0.9, 1.45)
    for idx, color in enumerate(colors):
        for row_idx in (0, 1, 2):
            cell = table[row_idx, idx]
            cell.get_text().set_color(color)
            _apply_text_outline(cell.get_text(), linewidth=0.3)
    for cell in table.get_celld().values():
        cell.set_linewidth(0.0)
        cell.set_edgecolor("none")
        cell.set_facecolor("none")
    ax.set_position([0.1, 0.52, 0.8, 0.4])
    fig.tight_layout()
    out_path = out_dir / f"{output_prefix}_filman_walls_composition_mass.png"
    _set_transparent(fig, ax)
    fig.savefig(out_path, dpi=plot_dpi, bbox_inches="tight", transparent=True)
    plt.close(fig)


def write_filman_walls_violin_plot(
    out_dir: Path,
    output_prefix: str,
    values_by_cat: dict[str, list[float]],
    have_plots: bool,
    scalar_name: str,
    font_scale: float,
    plot_dpi: int,
    plot_percentile_range: tuple[float, float] | None = None,
    log10_field_value: bool = False,
) -> None:
    if not have_plots:
        return
    try:
        import numpy as np
    except Exception:
        return

    data: list[np.ndarray] = []
    keys, labels, colors = _composition_lists()
    stats_by_cat = []
    for key in keys:
        vals = values_by_cat.get(key, [])
        arr = np.asarray(vals, dtype=float) if len(vals) else np.array([], dtype=float)
        if plot_percentile_range:
            arr = trim_by_percentile(arr, plot_percentile_range[0], plot_percentile_range[1])
        data.append(arr)
        stats_by_cat.append(aggregate_values_np(arr))

    import matplotlib.pyplot as plt

    sizes = _font_sizes(font_scale)
    fig, ax = plt.subplots(figsize=(9, 4.2))
    parts = ax.violinplot(data, showmeans=False, showmedians=True, showextrema=False)
    for body, color in zip(parts["bodies"], colors):
        body.set_facecolor(color)
        body.set_edgecolor("black")
        body.set_alpha(1.0)
    ax.set_xticks([])
    ax.set_xticklabels([])
    label = _scalar_label(scalar_name)
    y_label = label
    title_label = label
    ax.set_ylabel(y_label, fontsize=sizes["label"])
    mins = [s["min"] for s in stats_by_cat if s["count"]]
    maxs = [s["max"] for s in stats_by_cat if s["count"]]
    if mins and maxs:
        ax.set_ylim(min(mins), max(maxs))
    ax.set_title(f"{title_label} Distribution by Category", fontsize=sizes["title"])
    ax.tick_params(axis="y", labelsize=sizes["tick"])

    stats_rows = ["max", "q75", "mean", "median", "q25", "min", "sd", "skew", "count"]
    cell_text = [
        [_fmt_table_num(s["max"]) if s["count"] else "" for s in stats_by_cat],
        [_fmt_table_num(s["q75"]) if s["count"] else "" for s in stats_by_cat],
        [_fmt_table_num(s["mean"]) if s["count"] else "" for s in stats_by_cat],
        [_fmt_table_num(s["median"]) if s["count"] else "" for s in stats_by_cat],
        [_fmt_table_num(s["q25"]) if s["count"] else "" for s in stats_by_cat],
        [_fmt_table_num(s["min"]) if s["count"] else "" for s in stats_by_cat],
        [_fmt_table_num(s["std"]) if s["count"] else "" for s in stats_by_cat],
        [_fmt_table_num(s["skew"]) if s["count"] else "" for s in stats_by_cat],
        [_fmt_table_sci(s["count"]) if s["count"] else "" for s in stats_by_cat],
    ]
    col_widths = [1.0 / len(labels)] * len(labels)
    table = ax.table(
        cellText=cell_text,
        rowLabels=stats_rows,
        colLabels=labels,
        cellLoc="center",
        rowLoc="center",
        loc="bottom",
        bbox=[0.0, -0.8, 1.0, 0.7],
        colWidths=col_widths,
    )
    table.auto_set_font_size(False)
    table.set_fontsize(sizes["table"])
    for idx, color in enumerate(colors):
        for row_idx in range(len(stats_rows) + 1):
            cell = table[row_idx, idx]
            cell.get_text().set_color(color)
            _apply_text_outline(cell.get_text(), linewidth=0.3)
            cell.get_text().set_ha("center")
    for cell in table.get_celld().values():
        cell.get_text().set_ha("center")
    for cell in table.get_celld().values():
        cell.set_linewidth(0.0)
        cell.set_edgecolor("none")
        cell.set_facecolor("none")
    for spine_name, spine in ax.spines.items():
        spine.set_visible(spine_name == "left")
    ax.set_position([0.1, 0.48, 0.8, 0.42])
    fig.tight_layout()
    out_path = out_dir / f"{output_prefix}_filman_walls_{scalar_name}_violin.png"
    _set_transparent(fig, ax)
    fig.savefig(out_path, dpi=plot_dpi, bbox_inches="tight", transparent=True)
    plt.close(fig)


def run_polars_chunked(paths: list[Path], args: argparse.Namespace) -> None:
    import polars as pl
    violin_scalar = _plot_scalar_name(args.violin_scalar, args.log10_field_value)
    box_scalar = _plot_scalar_name(args.box_scalar, args.log10_field_value)

    if args.hist_bin_mode != "global":
        print("[warn] Chunked polars only supports global histogram bins; forcing --hist-bin-mode global.")
        args.hist_bin_mode = "global"

    scalars, include_filman, include_cluster = _infer_schema_from_headers(paths)
    if not scalars:
        raise SystemExit("[error] No scalar columns found in topology_points.csv inputs.")

    # Use all scalars for histogram bins so quantiles can be reconstructed.
    hist_scalars = scalars

    plow = phigh = None
    if args.hist_percentile_range:
        plow, phigh = args.hist_percentile_range
        plow = plow / 100.0
        phigh = phigh / 100.0
    qlo = 0.0 if plow is None else plow
    qhi = 1.0 if phigh is None else phigh

    chunks = _split_chunks(paths, max(1, args.polars_chunks))

    # Global min/max (approximate if percentile range requested).
    global_support: dict[str, tuple[float, float]] = {}
    for name in hist_scalars:
        low = float("inf")
        high = float("-inf")
        for chunk_paths in chunks:
            lazy = _prepare_polars_lazy(pl.concat([pl.scan_csv(str(p)) for p in chunk_paths]), include_filman, include_cluster)
            if plow is not None and phigh is not None:
                row = _collect(
                    lazy.select(
                        [
                            pl.col(name).quantile(plow, "linear").alias("min"),
                            pl.col(name).quantile(phigh, "linear").alias("max"),
                        ]
                    )
                ).to_dicts()[0]
            else:
                row = _collect(
                    lazy.select(
                        [pl.col(name).min().alias("min"), pl.col(name).max().alias("max")]
                    )
                ).to_dicts()[0]
            if row.get("min") is None or row.get("max") is None:
                continue
            low = min(low, float(row["min"]))
            high = max(high, float(row["max"]))
        if low != float("inf") and high != float("-inf") and high > low:
            global_support[name] = (low, high)

    if not global_support:
        raise SystemExit("[error] No valid scalar ranges found for histogram bins.")

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_keys, _, _ = _composition_lists()
    plot_scalars = {violin_scalar, box_scalar}
    plot_scalars = {s for s in plot_scalars if s in scalars}
    collect_raw_plots = bool(plot_scalars) and not args.no_plots
    plot_cache_dir = out_dir / f"{args.output_prefix}_plot_cache"
    plot_files: dict[str, dict[str, "io.BufferedWriter"]] = {}
    if collect_raw_plots:
        plot_cache_dir.mkdir(parents=True, exist_ok=True)
        for name in plot_scalars:
            plot_files[name] = {}
            for cat in plot_keys:
                bin_path = plot_cache_dir / f"{args.output_prefix}_{cat}_{name}.bin"
                plot_files[name][cat] = open(bin_path, "wb")

    # Initialize accumulators.
    total_counts: dict[str, int] = {}
    sums: dict[str, dict[str, float]] = {}
    sumsqs: dict[str, dict[str, float]] = {}
    mins: dict[str, dict[str, float]] = {}
    maxs: dict[str, dict[str, float]] = {}
    hist_counts: dict[str, dict[str, list[int]]] = {}
    hist_edges: dict[str, dict[str, list[float]]] = {}

    try:
        import numpy as np
    except Exception as exc:
        raise SystemExit(f"[error] numpy is required for histograms: {exc}")

    edges_by_scalar: dict[str, list[float]] = {}
    for name in hist_scalars:
        if name not in global_support:
            continue
        low, high = global_support[name]
        edges_by_scalar[name] = np.linspace(low, high, args.bins + 1, dtype=float).tolist()

    for chunk_paths in chunks:
        lazy = _prepare_polars_lazy(pl.concat([pl.scan_csv(str(p)) for p in chunk_paths]), include_filman, include_cluster)
        data = polars_category_frame(lazy, include_filman, include_cluster)

        agg_exprs = [pl.len().alias("count")]
        for name in scalars:
            agg_exprs.extend(
                [
                    pl.col(name).sum().alias(f"{name}_sum"),
                    (pl.col(name) * pl.col(name)).sum().alias(f"{name}_sumsq"),
                    pl.col(name).min().alias(f"{name}_min"),
                    pl.col(name).max().alias(f"{name}_max"),
                ]
            )
        stats_df = _collect(data.group_by("category").agg(agg_exprs))
        for row in stats_df.to_dicts():
            cat = row["category"]
            total_counts[cat] = total_counts.get(cat, 0) + int(row.get("count", 0) or 0)
            for name in scalars:
                sums.setdefault(cat, {})[name] = sums.get(cat, {}).get(name, 0.0) + float(row.get(f"{name}_sum", 0.0) or 0.0)
                sumsqs.setdefault(cat, {})[name] = sumsqs.get(cat, {}).get(name, 0.0) + float(row.get(f"{name}_sumsq", 0.0) or 0.0)
                mn = row.get(f"{name}_min", None)
                mx = row.get(f"{name}_max", None)
                if mn is not None:
                    mins.setdefault(cat, {})[name] = min(mins.get(cat, {}).get(name, float("inf")), float(mn))
                if mx is not None:
                    maxs.setdefault(cat, {})[name] = max(maxs.get(cat, {}).get(name, float("-inf")), float(mx))

        if collect_raw_plots:
            plot_data = data.filter(pl.col("category").is_in(plot_keys))
            for name in plot_scalars:
                grouped_vals = _collect(
                    plot_data.select(["category", pl.col(name)])
                    .group_by("category")
                    .agg(pl.col(name).alias("vals"))
                )
                for row in grouped_vals.to_dicts():
                    cat = row.get("category")
                    vals = row.get("vals", [])
                    if cat in plot_files.get(name, {}) and vals:
                        np.asarray(vals, dtype=float).tofile(plot_files[name][cat])

        for name in hist_scalars:
            if name not in edges_by_scalar:
                continue
            edges = edges_by_scalar[name]
            edges_arr = np.asarray(edges, dtype=float)

            def _bin_series(series: "pl.Series") -> "pl.Series":
                arr = series.to_numpy()
                idx = np.searchsorted(edges_arr, arr, side="right") - 1
                idx[(arr < edges_arr[0]) | (arr > edges_arr[-1])] = -1
                return pl.Series(idx, dtype=pl.Int64)

            bin_expr = pl.col(name).map_batches(_bin_series, return_dtype=pl.Int64).alias("bin")
            grouped = _collect(
                data.select(["category", bin_expr])
                .filter(pl.col("bin") >= 0)
                .group_by(["category", "bin"])
                .agg(pl.len().alias("count"))
            )
            for cat in total_counts:
                hist_counts.setdefault(cat, {}).setdefault(name, [0] * (len(edges) - 1))
                hist_edges.setdefault(cat, {}).setdefault(name, edges)
            for row in grouped.to_dicts():
                cat = row["category"]
                bin_idx = int(row["bin"])
                cnt = int(row["count"])
                if cat not in hist_counts:
                    hist_counts[cat] = {name: [0] * (len(edges) - 1)}
                    hist_edges[cat] = {name: edges}
                if name not in hist_counts[cat]:
                    hist_counts[cat][name] = [0] * (len(edges) - 1)
                    hist_edges[cat][name] = edges
                if 0 <= bin_idx < len(hist_counts[cat][name]):
                    hist_counts[cat][name][bin_idx] += cnt

    for name_map in plot_files.values():
        for handle in name_map.values():
            handle.close()

    exact_quantiles: dict[str, dict[str, tuple[float, float, float]]] = {}
    if collect_raw_plots and plot_scalars:
        try:
            import numpy as np
        except Exception:
            np = None
        if np is not None:
            for name in plot_scalars:
                raw_vals = plot_values_from_raw_files(
                    plot_cache_dir, args.output_prefix, name, plot_keys
                )
                for cat, arr in raw_vals.items():
                    if arr is None or len(arr) == 0:
                        continue
                    q25, median, q75 = np.quantile(arr, [0.25, 0.5, 0.75])
                    exact_quantiles.setdefault(cat, {})[name] = (
                        float(q25),
                        float(median),
                        float(q75),
                    )

    stats_path = out_dir / f"{args.output_prefix}_topology_stats.csv"
    fieldnames = _build_fieldnames(scalars)

    medians_map: dict[str, dict[str, float]] = {}
    box_stats_map: dict[str, dict[str, dict[str, float]]] = {}
    with open(stats_path, "w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        stats_rows: list[dict[str, object]] = []
        for cat in ALL_CATEGORIES:
            row: dict[str, object] = {"category": cat, "count": total_counts.get(cat, 0)}
            medians_map.setdefault(cat, {})
            for name in scalars:
                count = total_counts.get(cat, 0)
                sum_val = sums.get(cat, {}).get(name, 0.0)
                mean_val = sum_val / count if count else 0.0
                min_val = mins.get(cat, {}).get(name, 0.0)
                max_val = maxs.get(cat, {}).get(name, 0.0)
                sumsq = sumsqs.get(cat, {}).get(name, 0.0)
                var = (sumsq / count - mean_val * mean_val) if count else 0.0
                std_val = float(np.sqrt(var)) if count and var >= 0 else 0.0
                q25 = median = q75 = 0.0
                exact = exact_quantiles.get(cat, {}).get(name)
                if exact is not None:
                    q25, median, q75 = exact
                else:
                    edges = hist_edges.get(cat, {}).get(name)
                    counts = hist_counts.get(cat, {}).get(name)
                    if edges and counts:
                        q25 = _hist_quantile_from_counts(counts, edges, 0.25)
                        median = _hist_quantile_from_counts(counts, edges, 0.50)
                        q75 = _hist_quantile_from_counts(counts, edges, 0.75)
                row[f"{name}_sum"] = sum_val
                row[f"{name}_mean"] = mean_val
                row[f"{name}_min"] = min_val
                row[f"{name}_q25"] = q25
                row[f"{name}_median"] = median
                row[f"{name}_q75"] = q75
                row[f"{name}_max"] = max_val
                row[f"{name}_std"] = std_val
                medians_map[cat][name] = median
                box_stats_map.setdefault(name, {})[cat] = {
                    "count": float(count),
                    "mean": float(mean_val),
                    "min": float(min_val),
                    "q25": float(q25),
                    "median": float(median),
                    "q75": float(q75),
                    "max": float(max_val),
                    "std": float(std_val),
                }
            writer.writerow(row)
            stats_rows.append(row)
    print(f"[done] wrote {stats_path}")
    have_plots = _check_matplotlib(args.no_plots)
    cat_counts = total_counts
    _write_composition_and_transform_stats(out_dir, args.output_prefix, stats_rows, scalars, have_plots, args)
    _write_hist_plots(out_dir, args.output_prefix, hist_counts, hist_edges, medians_map, cat_counts, have_plots, args)

    if have_plots:
        plot_scalars = {violin_scalar, box_scalar}
        plot_scalars = {s for s in plot_scalars if s in scalars}
        if plot_scalars:
            plot_values: dict[str, dict[str, list[float]]] = {}
            focus_keys, _, _ = _composition_lists()
            if collect_raw_plots:
                for name in plot_scalars:
                    plot_values[name] = plot_values_from_raw_files(
                        plot_cache_dir,
                        args.output_prefix,
                        name,
                        focus_keys,
                    )
                if violin_scalar in plot_values:
                    write_filman_walls_violin_plot(
                        out_dir,
                        args.output_prefix,
                        plot_values[violin_scalar],
                        have_plots,
                        scalar_name=violin_scalar,
                        font_scale=args.plot_fontscale,
                        plot_dpi=args.plot_dpi,
                        plot_percentile_range=tuple(args.plot_percentile_range)
                        if args.plot_percentile_range
                        else None,
                        log10_field_value=args.log10_field_value,
                    )
                if box_scalar in plot_values:
                    exact_box_stats: dict[str, dict[str, float]] = {}
                    for key in focus_keys:
                        vals = plot_values[box_scalar].get(key, [])
                        if args.plot_percentile_range:
                            vals = trim_by_percentile(
                                vals, args.plot_percentile_range[0], args.plot_percentile_range[1]
                            )
                        exact_box_stats[key] = aggregate_values_np(vals)
                    write_filman_walls_box_plot(
                        out_dir,
                        args.output_prefix,
                        exact_box_stats,
                        have_plots,
                        scalar_name=box_scalar,
                        font_scale=args.plot_fontscale,
                        plot_dpi=args.plot_dpi,
                        log10_field_value=args.log10_field_value,
                    )
            else:
                for name in plot_scalars:
                    plot_values[name] = plot_values_from_hist(
                        hist_counts, hist_edges, name, args.plot_sample_size
                    )
                    have_focus = [k for k in focus_keys if plot_values[name].get(k)]
                    if len(have_focus) < len(focus_keys):
                        fallback = plot_values_from_hist_csvs(
                            out_dir,
                            args.output_prefix,
                            name,
                            focus_keys,
                            args.plot_sample_size,
                            seed=0,
                        )
                        for key, vals in fallback.items():
                            if vals:
                                plot_values[name][key] = vals
                if violin_scalar in plot_values:
                    write_filman_walls_violin_plot(
                        out_dir,
                        args.output_prefix,
                        plot_values[violin_scalar],
                        have_plots,
                        scalar_name=violin_scalar,
                        font_scale=args.plot_fontscale,
                        plot_dpi=args.plot_dpi,
                        plot_percentile_range=tuple(args.plot_percentile_range)
                        if args.plot_percentile_range
                        else None,
                        log10_field_value=args.log10_field_value,
                    )
                if box_scalar in scalars:
                    if args.plot_percentile_range and box_scalar in plot_values:
                        exact_box_stats: dict[str, dict[str, float]] = {}
                        for key in focus_keys:
                            vals = plot_values[box_scalar].get(key, [])
                            vals = trim_by_percentile(
                                vals, args.plot_percentile_range[0], args.plot_percentile_range[1]
                            )
                            exact_box_stats[key] = aggregate_values_np(vals)
                        write_filman_walls_box_plot(
                            out_dir,
                            args.output_prefix,
                            exact_box_stats,
                            have_plots,
                            scalar_name=box_scalar,
                            font_scale=args.plot_fontscale,
                            plot_dpi=args.plot_dpi,
                            log10_field_value=args.log10_field_value,
                        )
                    else:
                        write_filman_walls_box_plot(
                            out_dir,
                            args.output_prefix,
                            box_stats_map.get(box_scalar, {}),
                            have_plots,
                            scalar_name=box_scalar,
                            font_scale=args.plot_fontscale,
                            plot_dpi=args.plot_dpi,
                            log10_field_value=args.log10_field_value,
                        )
def write_filman_walls_box_plot(
    out_dir: Path,
    output_prefix: str,
    stats_by_cat: dict[str, dict[str, float]],
    have_plots: bool,
    scalar_name: str,
    font_scale: float,
    plot_dpi: int,
    log10_field_value: bool = False,
) -> None:
    if not have_plots:
        return

    import matplotlib.pyplot as plt

    sizes = _font_sizes(font_scale)
    keys, labels, colors = _composition_lists()

    stats_list: list[dict[str, float]] = []
    for key in keys:
        stats_list.append(stats_by_cat.get(key, {"count": 0.0}))

    bxp_stats = []
    for stats in stats_list:
        if not stats or stats.get("count", 0) == 0:
            bxp_stats.append(
                {"med": 0.0, "q1": 0.0, "q3": 0.0, "whislo": 0.0, "whishi": 0.0, "fliers": []}
            )
            continue
        q1 = stats["q25"]
        q3 = stats["q75"]
        iqr = q3 - q1
        whisker_lo = max(stats["min"], q1 - 1.5 * iqr)
        whisker_hi = min(stats["max"], q3 + 1.5 * iqr)
        bxp_stats.append(
            {
                "med": stats["median"],
                "q1": q1,
                "q3": q3,
                "whislo": whisker_lo,
                "whishi": whisker_hi,
                "fliers": [],
            }
        )

    fig, ax = plt.subplots(figsize=(9, 4.2))
    bp = ax.bxp(
        bxp_stats,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black"},
        boxprops={"edgecolor": "black"},
        whiskerprops={"color": "black"},
        capprops={"color": "black"},
    )

    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(1.0)
    ax.set_xticks([])
    ax.set_xticklabels([])
    # ax.set_ylabel(scalar_name)
    label = _scalar_label(scalar_name)
    y_label = label
    title_label = label
    ax.set_ylabel(y_label, fontsize=sizes["label"])
    mins = [b["whislo"] for b in bxp_stats]
    maxs = [b["whishi"] for b in bxp_stats]
    if mins and maxs:
        ymin = min(mins)
        ymax = max(maxs)
        pad = (ymax - ymin) * 0.05 if ymax > ymin else 0.0
        ax.set_ylim(ymin, ymax + pad)
    # ax.set_title(f"{scalar_name} distribution (boxplot)")
    ax.set_title(f"{title_label} Distribution by Category", fontsize=sizes["title"])
    ax.tick_params(axis="y", labelsize=sizes["tick"])

    stats_rows = ["max", "q75", "mean", "median", "q25", "min", "sd", "skew", "count"]
    def _skew_nonparam(stats: dict[str, float]) -> float | None:
        if not stats or not stats.get("count"):
            return None
        mean = stats.get("mean", 0.0)
        median = stats.get("median", 0.0)
        std = stats.get("std", 0.0)
        if std == 0:
            return None
        return (mean - median) / std

    cell_text = [
        [_fmt_table_num(s.get("max", 0.0)) if s.get("count") else "" for s in stats_list],
        [_fmt_table_num(s.get("q75", 0.0)) if s.get("count") else "" for s in stats_list],
        [_fmt_table_num(s.get("mean", 0.0)) if s.get("count") else "" for s in stats_list],
        [_fmt_table_num(s.get("median", 0.0)) if s.get("count") else "" for s in stats_list],
        [_fmt_table_num(s.get("q25", 0.0)) if s.get("count") else "" for s in stats_list],
        [_fmt_table_num(s.get("min", 0.0)) if s.get("count") else "" for s in stats_list],
        [_fmt_table_num(s.get("std", 0.0)) if s.get("count") else "" for s in stats_list],
        [
            _fmt_table_num(val) if val is not None else ""
            for val in (_skew_nonparam(s) for s in stats_list)
        ],
        [_fmt_table_sci(s.get("count", 0.0)) if s.get("count") else "" for s in stats_list],
    ]
    col_widths = [1.0 / len(labels)] * len(labels)
    table = ax.table(
        cellText=cell_text,
        rowLabels=stats_rows,
        colLabels=labels,
        cellLoc="center",
        rowLoc="center",
        loc="bottom",
        bbox=[0.0, -0.8, 1.0, 0.7],
        colWidths=col_widths,
    )
    table.auto_set_font_size(False)
    table.set_fontsize(sizes["table"])
    col_colors = _composition_color_map()
    for (row_idx, col_idx), cell in table.get_celld().items():
        if col_idx in col_colors:
            cell.get_text().set_color(col_colors[col_idx])
            _apply_text_outline(cell.get_text(), linewidth=0.3)
        if col_idx in col_colors or row_idx == 0:
            cell.get_text().set_ha("center")
    for cell in table.get_celld().values():
        cell.get_text().set_ha("center")
    for cell in table.get_celld().values():
        cell.set_linewidth(0.0)
        cell.set_edgecolor("none")
        cell.set_facecolor("none")
    for spine_name, spine in ax.spines.items():
        spine.set_visible(spine_name == "left")
    ax.set_position([0.1, 0.48, 0.8, 0.42])
    fig.tight_layout()
    out_path = out_dir / f"{output_prefix}_filman_walls_{scalar_name}_box.png"
    _set_transparent(fig, ax)
    fig.savefig(out_path, dpi=plot_dpi, bbox_inches="tight", transparent=True)
    plt.close(fig)


ALL_CATEGORIES = [
    # Base categories
    "walls",
    "walls_only",
    "filmans",
    "filmans_only",
    "clusters",
    "clusters_only",
    "unassigned",
    # Exclusions
    "walls_not_clusters",
    "walls_not_filmans",
    "filmans_not_clusters",
    "filmans_not_walls",
    "clusters_not_filmans",
    "clusters_not_walls",
    # Universe
    "unassigned_walls_filmans_clusters",
    # Exclusive pairwise intersections (sigma-algebra atoms)
    "shared_walls_filmans",
    "shared_walls_clusters",
    "shared_filmans_clusters",
    "shared_walls_filmans_clusters",
    # Sigma-algebra: inclusive intersections
    "walls_and_filmans",
    "walls_and_clusters",
    "filmans_and_clusters",
    # Pairwise unions
    "walls_or_filmans",
    "walls_or_clusters",
    "filmans_or_clusters",
    "all_structures",
    # Complements of generators
    "not_walls",
    "not_filmans",
    "not_clusters",
    # Complements of pairwise intersections
    "not_walls_and_filmans",
    "not_walls_and_clusters",
    "not_filmans_and_clusters",
    # Complement of triple intersection
    "not_walls_and_filmans_and_clusters",
    # Complements of pairwise unions
    "neither_walls_nor_filmans",
    "neither_walls_nor_clusters",
    "neither_filmans_nor_clusters",
]


COMPOSITION_COMPONENTS = [
    ("clusters", "Clusters", "#f2e661"),
    ("filmans_not_clusters", "Filaments", "#fca50a"),
    ("walls_only", "Walls", "#c73e4c"),
    ("unassigned", "Unassigned", "#5d126e"),
]


def main() -> None:
    args = parse_args()
    violin_scalar = _plot_scalar_name(args.violin_scalar, args.log10_field_value)
    box_scalar = _plot_scalar_name(args.box_scalar, args.log10_field_value)
    if args.inputs:
        paths = [Path(p) for p in args.inputs]
    else:
        root = Path(args.root)
        paths = sorted(root.glob(args.glob))
    if not paths:
        raise SystemExit("[error] No topology_points.csv files found.")

    if args.engine == "polars":
        try:
            import polars as pl  # noqa: F401
        except Exception as exc:
            print(f"[warn] Polars unavailable ({exc}); falling back to python engine.")
            args.engine = "python"

    if args.engine == "polars":
        if args.polars_chunks and args.polars_chunks > 1:
            run_polars_chunked(paths, args)
        else:
            run_polars(paths, args)
        return
    if args.engine == "stream":
        run_stream(paths, args)
        return

    scalars, values = read_inputs(paths, args.log10_field_value)
    if not scalars:
        raise SystemExit("[error] No scalar columns found in topology_points.csv inputs.")
    hist_scalars = _select_hist_scalars(scalars, args)

    cat_counts: dict[str, int] = {}
    for cat, cat_vals in values.items():
        count = 0
        for name in scalars:
            if name in cat_vals:
                count = len(cat_vals[name])
                break
        cat_counts[cat] = count

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    stats_path = out_dir / f"{args.output_prefix}_topology_stats.csv"

    fieldnames = _build_fieldnames(scalars)

    stats_rows: list[dict[str, object]] = []
    with open(stats_path, "w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        box_stats_map: dict[str, dict[str, dict[str, float]]] = {}
        for cat in ALL_CATEGORIES:
            row: dict[str, object] = {"category": cat}
            for name in scalars:
                stats = aggregate_values(values.get(cat, {}).get(name, []))
                row[f"{name}_sum"] = stats["sum"]
                row[f"{name}_mean"] = stats["mean"]
                row[f"{name}_min"] = stats["min"]
                row[f"{name}_q25"] = stats["q25"]
                row[f"{name}_median"] = stats["median"]
                row[f"{name}_q75"] = stats["q75"]
                row[f"{name}_max"] = stats["max"]
                row[f"{name}_std"] = stats["std"]
                row["count"] = stats["count"]
                box_stats_map.setdefault(name, {})[cat] = stats
            writer.writerow(row)
            stats_rows.append(row)
    print(f"[done] wrote {stats_path}")
    _write_transform_stats(out_dir, args.output_prefix, stats_rows, scalars)

    # Histograms
    try:
        import numpy as np
    except Exception as exc:
        raise SystemExit(f"[error] numpy is required for histograms: {exc}")

    have_plots = _check_matplotlib(args.no_plots)
    composition_counts = _composition_counts_from_rows(stats_rows)
    write_filman_walls_composition_plot(
        out_dir,
        args.output_prefix,
        composition_counts,
        have_plots,
        args.plot_fontscale,
        args.plot_dpi,
    )

    edges_map: dict[str, dict[str, list[float]]] = {}
    plow = phigh = None
    if args.hist_percentile_range:
        plow, phigh = args.hist_percentile_range
        plow = plow / 100.0
        phigh = phigh / 100.0
    qlo = 0.0 if plow is None else plow
    qhi = 1.0 if phigh is None else phigh

    global_support: dict[str, tuple[float, float]] = {}
    for name in hist_scalars:
        all_vals: list[float] = []
        for cat_vals in values.values():
            all_vals.extend(cat_vals.get(name, []))
        if not all_vals:
            continue
        vals_sorted = sorted(all_vals)
        if plow is not None and phigh is not None:
            gmin = _quantile(vals_sorted, plow)
            gmax = _quantile(vals_sorted, phigh)
        else:
            gmin = float(vals_sorted[0])
            gmax = float(vals_sorted[-1])
        if gmax > gmin:
            global_support[name] = (gmin, gmax)

    if args.hist_bin_mode == "global":
        for name in hist_scalars:
            if name not in global_support:
                continue
            gmin, gmax = global_support[name]
            if gmax <= gmin:
                continue
            edges = [gmin + (gmax - gmin) * i / args.bins for i in range(args.bins + 1)]
            for cat in values:
                edges_map.setdefault(cat, {})[name] = edges
    else:
        for cat, cat_vals in values.items():
            for name in hist_scalars:
                vals = cat_vals.get(name, [])
                if not vals or name not in global_support:
                    continue
                vals_sorted = sorted(vals)
                gmin, gmax = global_support[name]
                edges = _category_quantile_edges(vals_sorted, args.bins, qlo, qhi, gmin, gmax)
                if edges:
                    edges_map.setdefault(cat, {})[name] = edges

    hist_counts: dict[str, dict[str, list[int]]] = {}
    hist_edges: dict[str, dict[str, list[float]]] = {}
    for cat, cat_vals in values.items():
        for name in hist_scalars:
            vals = cat_vals.get(name, [])
            if not vals:
                continue
            if name not in edges_map.get(cat, {}):
                continue
            edges = edges_map[cat][name]
            arr = np.asarray(vals, dtype=float)
            counts, edges = np.histogram(arr, bins=edges)
            hist_counts.setdefault(cat, {})[name] = counts.tolist()
            hist_edges.setdefault(cat, {})[name] = edges.tolist()
            csv_path = out_dir / f"{args.output_prefix}_{cat}_{name}_hist.csv"
            write_histogram_csv(csv_path, edges.tolist(), counts.tolist(), len(vals))
            if have_plots:
                import matplotlib.pyplot as plt
                centers = [(edges[i] + edges[i + 1]) / 2 for i in range(len(edges) - 1)]
                widths = [edges[i + 1] - edges[i] for i in range(len(edges) - 1)]
                fig = plt.figure()
                plt.bar(centers, counts, width=widths, align="center")
                plt.axvline(float(np.median(arr)), color="red", linestyle="--", linewidth=1)
                plt.xlabel(_scalar_label(name))
                plt.ylabel("count")
                plt.title(f"{cat}: {_scalar_label(name)} (n={len(vals):.2e})")
                fig_path = out_dir / f"{args.output_prefix}_{cat}_{name}_hist.png"
                fig.savefig(fig_path, dpi=args.plot_dpi, bbox_inches="tight")
                plt.close(fig)
    if have_plots:
        plot_scalars = {violin_scalar, box_scalar}
        plot_scalars = {s for s in plot_scalars if s in scalars}
        if not args.plots_from_raw and hist_counts and hist_edges:
            plot_values = {}
            for name in plot_scalars:
                if name in hist_scalars:
                    plot_values[name] = plot_values_from_hist(
                        hist_counts, hist_edges, name, args.plot_sample_size
                    )
                else:
                    plot_values[name] = {}
                    for cat in values:
                        plot_values[name][cat] = values[cat].get(name, [])
        else:
            plot_values = {name: {} for name in plot_scalars}
            for name in plot_scalars:
                for cat in values:
                    plot_values[name][cat] = values[cat].get(name, [])
        if violin_scalar in plot_values:
            write_filman_walls_violin_plot(
                out_dir,
                args.output_prefix,
                plot_values[violin_scalar],
                have_plots,
                scalar_name=violin_scalar,
                font_scale=args.plot_fontscale,
                plot_dpi=args.plot_dpi,
                plot_percentile_range=tuple(args.plot_percentile_range)
                if args.plot_percentile_range
                else None,
                log10_field_value=args.log10_field_value,
            )
        if box_scalar in plot_values:
            if args.plot_percentile_range:
                exact_box_stats: dict[str, dict[str, float]] = {}
                keys, _, _ = _composition_lists()
                for key in keys:
                    vals = plot_values[box_scalar].get(key, [])
                    vals = trim_by_percentile(
                        vals, args.plot_percentile_range[0], args.plot_percentile_range[1]
                    )
                    exact_box_stats[key] = aggregate_values_np(vals)
                write_filman_walls_box_plot(
                    out_dir,
                    args.output_prefix,
                    exact_box_stats,
                    have_plots,
                    scalar_name=box_scalar,
                    font_scale=args.plot_fontscale,
                    plot_dpi=args.plot_dpi,
                    log10_field_value=args.log10_field_value,
                )
            else:
                write_filman_walls_box_plot(
                    out_dir,
                    args.output_prefix,
                    box_stats_map.get(box_scalar, {}),
                    have_plots,
                    scalar_name=box_scalar,
                    font_scale=args.plot_fontscale,
                    plot_dpi=args.plot_dpi,
                    log10_field_value=args.log10_field_value,
                )



if __name__ == "__main__":
    main()
