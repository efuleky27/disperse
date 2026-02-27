#!/usr/bin/env python3
"""
Aggregate per-point topology CSVs across crops into summary stats and histograms.

Reads *_topology_points.csv files (from ndtopo_stats.py --per-point-csv) and
computes the same category stats as topology_stats.csv, plus histogram outputs.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


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
        "--hist-percentile-range",
        type=float,
        nargs=2,
        metavar=("PLOW", "PHIGH"),
        help="Percentile range for histogram bounds (e.g., 1 99).",
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
        default=int(10e+8),
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
    return p.parse_args()


def _quantile(vals_sorted: List[float], q: float) -> float:
    n = len(vals_sorted)
    if n == 1:
        return float(vals_sorted[0])
    idx = q * (n - 1)
    lo = int(idx)
    hi = min(lo + 1, n - 1)
    frac = idx - lo
    return float(vals_sorted[lo] * (1 - frac) + vals_sorted[hi] * frac)




def _select_hist_scalars(all_scalars: List[str], args: argparse.Namespace) -> List[str]:
    if args.hist_scalars:
        if len(args.hist_scalars) == 1 and args.hist_scalars[0].lower() == "all":
            return list(all_scalars)
        return [name for name in args.hist_scalars if name in all_scalars]
    defaults = ["field_value", "log_field_value"]
    return [name for name in defaults if name in all_scalars]


def _font_sizes(scale: float) -> Dict[str, float]:
    base = {"title": 12.0, "label": 11.0, "table": 10.0, "tick": 10.0}
    return {key: val * scale for key, val in base.items()}


def _apply_text_outline(text_obj, linewidth: float = 1.0) -> None:
    try:
        import matplotlib.patheffects as pe
    except Exception:
        return
    text_obj.set_path_effects([pe.withStroke(linewidth=linewidth, foreground="black")])


def _component_lists() -> Tuple[List[str], List[str], List[str]]:
    keys: List[str] = []
    labels: List[str] = []
    colors: List[str] = []
    for key, label, color in FILMAN_WALLS_COMPONENTS:
        keys.append(key)
        labels.append(label)
        colors.append(color)
    return keys, labels, colors


def _component_color_map() -> Dict[int, str]:
    return {idx: color for idx, (_, _, color) in enumerate(FILMAN_WALLS_COMPONENTS)}


def _set_transparent(fig, ax) -> None:
    try:
        fig.patch.set_alpha(0.0)
    except Exception:
        pass
    try:
        ax.patch.set_alpha(0.0)
    except Exception:
        pass


def aggregate_values(vals: List[float]) -> Dict[str, float]:
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
        }
    vals_sorted = sorted(vals)
    n = len(vals_sorted)
    total = float(sum(vals_sorted))
    mean_val = total / n
    return {
        "count": n,
        "sum": total,
        "mean": mean_val,
        "min": float(vals_sorted[0]),
        "q25": _quantile(vals_sorted, 0.25),
        "median": _quantile(vals_sorted, 0.5),
        "q75": _quantile(vals_sorted, 0.75),
        "max": float(vals_sorted[-1]),
        "std": float((sum((v - mean_val) ** 2 for v in vals_sorted) / n) ** 0.5),
    }


def _normalize_edges(edges: List[float]) -> List[float]:
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
    vals_sorted: List[float], bins: int, plow: float, phigh: float
) -> List[float]:
    if not vals_sorted:
        return []
    qs = [plow + (phigh - plow) * i / bins for i in range(bins + 1)]
    edges = [_quantile(vals_sorted, q) for q in qs]
    return _normalize_edges(edges)


def _category_quantile_edges(
    vals_sorted: List[float],
    bins: int,
    plow: float,
    phigh: float,
    global_lo: float,
    global_hi: float,
) -> List[float]:
    if not vals_sorted:
        return []
    edges = _quantile_edges_from_sorted(vals_sorted, bins, plow, phigh)
    if not edges:
        return []
    edges[0] = float(global_lo)
    edges[-1] = float(global_hi)
    return _normalize_edges(edges)


def _bin_idx(edges: List[float], val: float) -> int:
    import bisect

    if not edges or val < edges[0] or val > edges[-1]:
        return -1
    idx = bisect.bisect_right(edges, val) - 1
    if idx < 0:
        return -1
    if idx >= len(edges) - 1:
        return len(edges) - 2
    return idx


def _sample_from_hist(
    edges: List[float],
    counts: List[int],
    sample_size: int,
    rng: "np.random.Generator",
) -> List[float]:
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
    hist_counts: Dict[str, Dict[str, List[int]]],
    hist_edges: Dict[str, Dict[str, List[float]]],
    scalar_name: str,
    sample_size: int,
    seed: int = 0,
) -> Dict[str, List[float]]:
    try:
        import numpy as np
    except Exception:
        return {}
    rng = np.random.default_rng(seed)
    out: Dict[str, List[float]] = {}
    for cat, name_map in hist_counts.items():
        counts = name_map.get(scalar_name)
        edges = hist_edges.get(cat, {}).get(scalar_name)
        if not counts or not edges:
            continue
        vals = _sample_from_hist(edges, counts, sample_size, rng)
        if vals:
            out[cat] = vals
    return out


def read_inputs(paths: Iterable[Path]) -> Tuple[List[str], Dict[str, Dict[str, List[float]]]]:
    scalars: List[str] = []
    hist_scalars: List[str] = []
    values: Dict[str, Dict[str, List[float]]] = {}
    for path in paths:
        with open(path, newline="", encoding="ascii") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames:
                continue
            if not scalars:
                scalar_candidates = [
                    name
                    for name in reader.fieldnames
                    if name not in ("delaunay_id", "is_wall", "is_filament", "is_filament_manifold")
                ]
                scalars = scalar_candidates
            for row in reader:
                is_wall = row.get("is_wall", "0").strip() == "1"
                is_fil = row.get("is_filament", "0").strip() == "1"
                is_filman = row.get("is_filament_manifold", "0").strip() == "1"

                categories: List[str] = []
                categories.append("walls") if is_wall else None
                categories.append("filaments") if is_fil else None
                if is_wall and is_fil:
                    categories.append("shared")
                if is_wall and not is_fil:
                    categories.append("wall_only")
                if is_fil and not is_wall:
                    categories.append("filament_only")
                if not is_wall and not is_fil:
                    categories.append("unassigned")
                if "is_filament_manifold" in row:
                    if is_filman:
                        categories.append("filament_manifolds")
                    if is_filman and not is_wall and not is_fil:
                        categories.append("filament_manifolds_only")
                    if is_filman and is_wall:
                        categories.append("shared_filament_manifolds_walls")
                    if is_filman and is_fil:
                        categories.append("shared_filaments_filament_manifolds")
                    if is_wall and not is_filman:
                        categories.append("walls_not_filament_manifolds")
                    if not is_wall and not is_filman:
                        categories.append("unassigned_walls_filament_manifolds")

                for cat in categories:
                    values.setdefault(cat, {})
                    for name in scalars:
                        if name not in row:
                            continue
                        try:
                            val = float(row[name])
                        except ValueError:
                            continue
                        values[cat].setdefault(name, []).append(val)
    return scalars, values


def polars_category_frame(lazy, include_filman: bool) -> "pl.LazyFrame":
    import polars as pl

    is_wall = pl.col("is_wall") == 1
    is_fil = pl.col("is_filament") == 1
    is_filman = pl.col("is_filament_manifold") == 1 if include_filman else pl.lit(False)

    cat_exprs = [
        pl.when(is_wall).then(pl.lit("walls")),
        pl.when(is_fil).then(pl.lit("filaments")),
        pl.when(is_wall & is_fil).then(pl.lit("shared")),
        pl.when(is_wall & ~is_fil).then(pl.lit("wall_only")),
        pl.when(is_fil & ~is_wall).then(pl.lit("filament_only")),
        pl.when(~is_wall & ~is_fil).then(pl.lit("unassigned")),
    ]
    if include_filman:
        cat_exprs.extend(
            [
                pl.when(is_filman).then(pl.lit("filament_manifolds")),
                pl.when(is_filman & ~is_wall & ~is_fil).then(pl.lit("filament_manifolds_only")),
                pl.when(is_filman & is_wall).then(pl.lit("shared_filament_manifolds_walls")),
                pl.when(is_filman & is_fil).then(pl.lit("shared_filaments_filament_manifolds")),
                pl.when(is_wall & ~is_filman).then(pl.lit("walls_not_filament_manifolds")),
                pl.when(~is_wall & ~is_filman).then(pl.lit("unassigned_walls_filament_manifolds")),
            ]
        )
    category_col = pl.concat_list(cat_exprs).list.drop_nulls().alias("category")
    return lazy.with_columns(category_col).explode("category")


def run_polars(paths: List[Path], args: argparse.Namespace) -> None:
    import polars as pl

    def _collect(lazy_frame: "pl.LazyFrame") -> "pl.DataFrame":
        try:
            return lazy_frame.collect(engine="streaming")
        except TypeError:
            return lazy_frame.collect(streaming=True)

    lazy_frames = [pl.scan_csv(str(p)) for p in paths]
    lazy = pl.concat(lazy_frames)
    schema = lazy.collect_schema()
    cols = set(schema.names())
    include_filman = "is_filament_manifold" in cols
    if not include_filman:
        lazy = lazy.with_columns(pl.lit(0).alias("is_filament_manifold"))
    scalars = [
        name
        for name in schema.names()
        if name not in ("delaunay_id", "is_wall", "is_filament", "is_filament_manifold")
    ]
    if not scalars:
        raise SystemExit("[error] No scalar columns found in topology_points.csv inputs.")

    data = polars_category_frame(lazy, include_filman)
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
    stats_df = stats_df.select(fieldnames)
    stats_df.write_csv(stats_path)
    print(f"[done] wrote {stats_path}")
    cat_counts = {row["category"]: int(row["count"]) for row in stats_df.to_dicts()}
    medians_map: Dict[str, Dict[str, float]] = {}
    box_stats_map: Dict[str, Dict[str, Dict[str, float]]] = {}
    for row in stats_df.to_dicts():
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
            }


    if args.no_plots:
        have_plots = False
    else:
        try:
            import matplotlib.pyplot as plt  # noqa: F401

            have_plots = True
        except Exception:
            have_plots = False
    write_filman_walls_composition_plot(
        out_dir,
        args.output_prefix,
        cat_counts,
        have_plots,
        args.plot_fontscale,
        args.plot_dpi,
    )

    hist_scalars = _select_hist_scalars(scalars, args)
    plow = phigh = None
    if args.hist_percentile_range:
        plow, phigh = args.hist_percentile_range
        plow = plow / 100.0
        phigh = phigh / 100.0
    qlo = 0.0 if plow is None else plow
    qhi = 1.0 if phigh is None else phigh

    hist_counts: Dict[str, Dict[str, List[int]]] = {}
    hist_edges: Dict[str, Dict[str, List[float]]] = {}
    plot_values_all: Dict[str, Dict[str, List[float]]] = {}

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

            for cat, counts in hist_counts.items():
                if name not in counts:
                    continue
                csv_path = out_dir / f"{args.output_prefix}_{cat}_{name}_hist.csv"
                write_histogram_csv(csv_path, edges, counts[name], cat_counts.get(cat, 0))
                if have_plots:
                    import matplotlib.pyplot as plt

                    fig = plt.figure()
                    centers = [(edges[i] + edges[i + 1]) / 2 for i in range(len(edges) - 1)]
                    widths = [edges[i + 1] - edges[i] for i in range(len(edges) - 1)]
                    plt.bar(centers, counts[name], width=widths, align="center")
                    median_val = medians_map.get(cat, {}).get(name)
                    if median_val is not None:
                        plt.axvline(median_val, color="red", linestyle="--", linewidth=1)
                    plt.xlabel(name)
                    plt.ylabel("count")
                    plt.title(f"{cat}: {name} (n={cat_counts.get(cat, 0):.2e})")
                    fig_path = out_dir / f"{args.output_prefix}_{cat}_{name}_hist.png"
                    fig.savefig(fig_path, dpi=args.plot_dpi, bbox_inches="tight")
                    plt.close(fig)
    else:
        # Per-category quantile bins with global support per scalar.
        global_support: Dict[str, Tuple[float, float]] = {}
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

        plot_values_all = collect_plot_values(paths, scalars)
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
                csv_path = out_dir / f"{args.output_prefix}_{cat}_{name}_hist.csv"
                write_histogram_csv(csv_path, edges, counts, cat_counts.get(cat, 0))
                if have_plots:
                    import matplotlib.pyplot as plt

                    fig = plt.figure()
                    centers = [(edges[i] + edges[i + 1]) / 2 for i in range(len(edges) - 1)]
                    widths = [edges[i + 1] - edges[i] for i in range(len(edges) - 1)]
                    plt.bar(centers, counts, width=widths, align="center")
                    median_val = medians_map.get(cat, {}).get(name)
                    if median_val is not None:
                        plt.axvline(median_val, color="red", linestyle="--", linewidth=1)
                    plt.xlabel(name)
                    plt.ylabel("count")
                    plt.title(f"{cat}: {name} (n={cat_counts.get(cat, 0):.2e})")
                    fig_path = out_dir / f"{args.output_prefix}_{cat}_{name}_hist.png"
                    fig.savefig(fig_path, dpi=args.plot_dpi, bbox_inches="tight")
                    plt.close(fig)

    if have_plots:
        plot_scalars = {args.violin_scalar, args.box_scalar}
        plot_scalars = {s for s in plot_scalars if s in scalars}
        if plot_scalars:
            plot_values: Dict[str, Dict[str, List[float]]] = {}
            if not args.plots_from_raw and hist_counts and hist_edges:
                for name in plot_scalars:
                    if name in hist_scalars:
                        plot_values[name] = plot_values_from_hist(
                            hist_counts, hist_edges, name, args.plot_sample_size
                        )
                    else:
                        if not plot_values_all:
                            plot_values_all = collect_plot_values(paths, list(plot_scalars))
                        plot_values[name] = plot_values_all.get(name, {})
            else:
                if not plot_values_all:
                    plot_values_all = collect_plot_values(paths, list(plot_scalars))
                plot_values = {name: plot_values_all.get(name, {}) for name in plot_scalars}
            if args.violin_scalar in plot_values:
                write_filman_walls_violin_plot(
                    out_dir,
                    args.output_prefix,
                    plot_values[args.violin_scalar],
                    have_plots,
                    scalar_name=args.violin_scalar,
                    font_scale=args.plot_fontscale,
                    plot_dpi=args.plot_dpi,
                )
            if args.box_scalar in scalars:
                write_filman_walls_box_plot(
                    out_dir,
                    args.output_prefix,
                    box_stats_map.get(args.box_scalar, {}),
                    have_plots,
                    scalar_name=args.box_scalar,
                    font_scale=args.plot_fontscale,
                    plot_dpi=args.plot_dpi,
                )


def _categories(is_wall: bool, is_fil: bool, is_filman: bool, include_filman: bool) -> List[str]:
    cats: List[str] = []
    if is_wall:
        cats.append("walls")
    if is_fil:
        cats.append("filaments")
    if is_wall and not is_fil:
        cats.append("walls_not_filaments")
    if is_fil and not is_wall:
        cats.append("filaments_not_walls")
    if is_wall and is_fil:
        cats.append("shared_walls_filaments")
    if include_filman:
        if is_filman:
            cats.append("filament_manifolds")
        if is_filman and not is_wall and not is_fil:
            cats.append("filament_manifolds_only")
        if is_filman and is_wall:
            cats.append("shared_filament_manifolds_walls")
        if is_fil and is_filman:
            cats.append("shared_filaments_filament_manifolds")
        if is_wall and not is_filman:
            cats.append("walls_not_filament_manifolds")
        if not is_wall and not is_filman:
            cats.append("unassigned_walls_filament_manifolds")
    if not is_wall and not is_fil and (not include_filman or not is_filman):
        cats.append("unassigned")
    return cats


def collect_plot_values(
    paths: Iterable[Path], scalars_needed: List[str]
) -> Dict[str, Dict[str, List[float]]]:
    values: Dict[str, Dict[str, List[float]]] = {name: {} for name in scalars_needed}
    scalars_set = set(scalars_needed)
    for path in paths:
        with open(path, newline="", encoding="ascii") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames:
                continue
            include_filman = "is_filament_manifold" in reader.fieldnames
            for row in reader:
                is_wall = row.get("is_wall", "0").strip() == "1"
                is_fil = row.get("is_filament", "0").strip() == "1"
                is_filman = (
                    include_filman and row.get("is_filament_manifold", "0").strip() == "1"
                )
                cats = _categories(is_wall, is_fil, is_filman, include_filman)
                for name in scalars_set:
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


def run_stream(paths: List[Path], args: argparse.Namespace) -> None:
    import math
    import random

    sample_size = 20000
    scalars: List[str] = []
    include_filman = False
    stats: Dict[str, Dict[str, Dict[str, float]]] = {}
    category_count: Dict[str, int] = {}
    samples: Dict[str, Dict[str, List[float]]] = {}
    global_samples: Dict[str, List[float]] = {}
    global_counts: Dict[str, int] = {}
    def _stat_store(cat: str, name: str) -> Dict[str, float]:
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

    def _sample_store(cat: str, name: str) -> List[float]:
        return samples.setdefault(cat, {}).setdefault(name, [])

    def _global_sample_store(name: str) -> List[float]:
        return global_samples.setdefault(name, [])

    # pass 1: streaming stats + samples
    for path in paths:
        with open(path, newline="", encoding="ascii") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames:
                continue
            if not scalars:
                scalars = [
                    name
                    for name in reader.fieldnames
                    if name not in ("delaunay_id", "is_wall", "is_filament", "is_filament_manifold")
                ]
                hist_scalars = _select_hist_scalars(scalars, args)
                include_filman = "is_filament_manifold" in reader.fieldnames
            for row in reader:
                is_wall = row.get("is_wall", "0").strip() == "1"
                is_fil = row.get("is_filament", "0").strip() == "1"
                is_filman = include_filman and row.get("is_filament_manifold", "0").strip() == "1"
                cats = _categories(is_wall, is_fil, is_filman, include_filman)
                for cat in cats:
                    category_count[cat] = category_count.get(cat, 0) + 1
                for name in scalars:
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

    global_support: Dict[str, Tuple[float, float]] = {}
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


    edges_by_cat: Dict[str, Dict[str, List[float]]] = {}
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
    counts: Dict[str, Dict[str, List[int]]] = {}
    for path in paths:
        with open(path, newline="", encoding="ascii") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames:
                continue
            for row in reader:
                is_wall = row.get("is_wall", "0").strip() == "1"
                is_fil = row.get("is_filament", "0").strip() == "1"
                is_filman = include_filman and row.get("is_filament_manifold", "0").strip() == "1"
                cats = _categories(is_wall, is_fil, is_filman, include_filman)
                for name in scalars:
                    try:
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
                        import bisect

                        bin_idx = bisect.bisect_right(edges, val) - 1
                        if bin_idx < 0 or bin_idx >= len(edges) - 1:
                            continue
                        counts.setdefault(cat, {}).setdefault(name, [0] * (len(edges) - 1))
                        counts[cat][name][bin_idx] += 1


    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    stats_path = out_dir / f"{args.output_prefix}_topology_stats.csv"
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
    with open(stats_path, "w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        median_map: Dict[str, Dict[str, float]] = {}
        box_stats_map: Dict[str, Dict[str, Dict[str, float]]] = {}
        for cat in sorted(stats.keys()):
            row: Dict[str, object] = {"category": cat, "count": category_count.get(cat, 0)}
            for name in scalars:
                store = stats[cat].get(name)
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
                median_map.setdefault(cat, {})[name] = qvals["median"]
                box_stats_map.setdefault(name, {})[cat] = {
                    "count": float(store["count"]),
                    "mean": float(store["mean"]),
                    "min": float(store["min"]),
                    "q25": float(qvals["q25"]),
                    "median": float(qvals["median"]),
                    "q75": float(qvals["q75"]),
                    "max": float(store["max"]),
                }
            writer.writerow(row)
    print(f"[done] wrote {stats_path}")

    if args.no_plots:
        have_plots = False
    else:
        try:
            import matplotlib.pyplot as plt  # noqa: F401

            have_plots = True
        except Exception:
            have_plots = False
    write_filman_walls_composition_plot(
        out_dir,
        args.output_prefix,
        cat_counts,
        have_plots,
        args.plot_fontscale,
        args.plot_dpi,
    )

    for cat, name_map in counts.items():
        for name, hist in name_map.items():
            edges = edges_by_cat.get(cat, {}).get(name)
            if not edges:
                continue
            csv_path = out_dir / f"{args.output_prefix}_{cat}_{name}_hist.csv"
            write_histogram_csv(csv_path, edges, hist, cat_counts.get(cat, 0))
            if have_plots:
                import matplotlib.pyplot as plt

                centers = [(edges[i] + edges[i + 1]) / 2 for i in range(len(edges) - 1)]
                widths = [edges[i + 1] - edges[i] for i in range(len(edges) - 1)]
                fig = plt.figure()
                plt.bar(centers, hist, width=widths, align="center")
                median_val = median_map.get(cat, {}).get(name)
                if median_val is not None:
                    plt.axvline(median_val, color="red", linestyle="--", linewidth=1)
                plt.xlabel(name)
                plt.ylabel("count")
                plt.title(f"{cat}: {name} (n={cat_counts.get(cat, 0):.2e})")
                fig_path = out_dir / f"{args.output_prefix}_{cat}_{name}_hist.png"
                fig.savefig(fig_path, dpi=args.plot_dpi, bbox_inches="tight")
                plt.close(fig)
    if have_plots:
        plot_scalars = {args.violin_scalar, args.box_scalar}
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
                        plot_values[name] = collect_plot_values(paths, [name]).get(name, {})
            else:
                plot_values = collect_plot_values(paths, sorted(plot_scalars))
            if args.violin_scalar in plot_values:
                write_filman_walls_violin_plot(
                    out_dir,
                    args.output_prefix,
                    plot_values[args.violin_scalar],
                    have_plots,
                    scalar_name=args.violin_scalar,
                    font_scale=args.plot_fontscale,
                    plot_dpi=args.plot_dpi,
                )
            if args.box_scalar in scalars:
                write_filman_walls_box_plot(
                    out_dir,
                    args.output_prefix,
                    box_stats_map.get(args.box_scalar, {}),
                    have_plots,
                    scalar_name=args.box_scalar,
                    font_scale=args.plot_fontscale,
                    plot_dpi=args.plot_dpi,
                )


def write_histogram_csv(
    path: Path,
    edges: List[float],
    counts: List[int],
    total_count: int,
) -> None:
    with open(path, "w", newline="", encoding="ascii") as handle:
        writer = csv.writer(handle)
        writer.writerow(["bin_left", "bin_right", "count", "total_count"])
        for i, count in enumerate(counts):
            writer.writerow([edges[i], edges[i + 1], count, total_count])


def _hist_quantile_from_counts(
    counts_list: List[int], edges_list: List[float], q: float
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


def _hist_stats_from_counts(counts_list: List[int], edges_list: List[float]) -> Dict[str, float]:
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
    cat_counts: Dict[str, int],
    have_plots: bool,
    font_scale: float,
    plot_dpi: int,
) -> None:
    if not have_plots:
        return
    categories = FILMAN_WALLS_COMPONENTS
    keys, labels, colors = _component_lists()
    total = sum(cat_counts.get(key, 0) for key in keys)
    if total <= 0:
        return
    percents = [cat_counts.get(key, 0) / total * 100.0 for key in keys]

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
    table_bbox = [0.13, -0.7, 0.85, 0.5]
    table = ax.table(
        cellText=[
            [f"{c:.2e}" for c in counts],
            [f"{p:.1f}%" for p in percents],
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
    table.scale(0.9, 1.2)
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


def write_filman_walls_violin_plot(
    out_dir: Path,
    output_prefix: str,
    values_by_cat: Dict[str, List[float]],
    have_plots: bool,
    scalar_name: str,
    font_scale: float,
    plot_dpi: int,
) -> None:
    if not have_plots:
        return
    try:
        import numpy as np
    except Exception:
        return

    data: List[np.ndarray] = []
    keys, labels, colors = _component_lists()
    stats_by_cat = []
    for key in keys:
        vals = values_by_cat.get(key, [])
        data.append(np.asarray(vals, dtype=float) if vals else np.array([], dtype=float))
        stats_by_cat.append(aggregate_values(vals))

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
    # ax.set_ylabel(scalar_name)
    ax.set_ylabel(
        "Log-density" if scalar_name == "log_field_value" else scalar_name,
        fontsize=sizes["label"],
    )
    mins = [s["min"] for s in stats_by_cat if s["count"]]
    maxs = [s["max"] for s in stats_by_cat if s["count"]]
    if mins and maxs:
        ax.set_ylim(min(mins), max(maxs))
    # ax.set_title(f"{scalar_name} distribution by category")
    ax.set_title(
        "Log-density distribution by category"
        if scalar_name == "log_field_value"
        else f"{scalar_name} distribution by category",
        fontsize=sizes["title"],
    )
    ax.tick_params(axis="y", labelsize=sizes["tick"])

    # Summary table under the violins
    stats_rows = ["count", "mean", "min", "q25", "median", "q75", "max"]
    cell_text = [
        [f"{s['count']:.2e}" if s["count"] else "" for s in stats_by_cat],
        [f"{s['mean']:.3f}" if s["count"] else "" for s in stats_by_cat],
        [f"{s['min']:.3f}" if s["count"] else "" for s in stats_by_cat],
        [f"{s['q25']:.3f}" if s["count"] else "" for s in stats_by_cat],
        [f"{s['median']:.3f}" if s["count"] else "" for s in stats_by_cat],
        [f"{s['q75']:.3f}" if s["count"] else "" for s in stats_by_cat],
        [f"{s['max']:.3f}" if s["count"] else "" for s in stats_by_cat],
    ]
    table = ax.table(
        cellText=cell_text,
        rowLabels=stats_rows,
        colLabels=labels,
        cellLoc="center",
        rowLoc="center",
        loc="bottom",
        bbox=[0.0, -0.6, 1.0, 0.5],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(sizes["table"])
    for idx, color in enumerate(colors):
        for row_idx in range(len(stats_rows) + 1):
            cell = table[row_idx, idx]
            cell.get_text().set_color(color)
            _apply_text_outline(cell.get_text(), linewidth=0.3)
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


def write_filman_walls_box_plot(
    out_dir: Path,
    output_prefix: str,
    stats_by_cat: Dict[str, Dict[str, float]],
    have_plots: bool,
    scalar_name: str,
    font_scale: float,
    plot_dpi: int,
) -> None:
    if not have_plots:
        return

    import matplotlib.pyplot as plt

    sizes = _font_sizes(font_scale)
    keys, labels, colors = _component_lists()

    stats_list: List[Dict[str, float]] = []
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
    ax.set_ylabel(
        "Density" if scalar_name == "field_value" else scalar_name,
        fontsize=sizes["label"],
    )
    mins = [b["whislo"] for b in bxp_stats]
    maxs = [b["whishi"] for b in bxp_stats]
    if mins and maxs:
        ymin = min(mins)
        ymax = max(maxs)
        pad = (ymax - ymin) * 0.05 if ymax > ymin else 0.0
        ax.set_ylim(ymin, ymax + pad)
    # ax.set_title(f"{scalar_name} distribution (boxplot)")
    ax.set_title(
        "Density distribution by category"
        if scalar_name == "field_value"
        else f"{scalar_name} distribution by category",
        fontsize=sizes["title"],
    )
    ax.tick_params(axis="y", labelsize=sizes["tick"])

    stats_rows = ["count", "mean", "min", "q25", "median", "q75", "max", "skew"]
    def _bowley(stats: Dict[str, float]) -> Optional[float]:
        if not stats or not stats.get("count"):
            return None
        q1 = stats.get("q25", 0.0)
        q2 = stats.get("median", 0.0)
        q3 = stats.get("q75", 0.0)
        denom = q3 - q1
        if denom == 0:
            return None
        return (q3 + q1 - 2.0 * q2) / denom

    cell_text = [
        [f"{s.get('count', 0.0):.2e}" if s.get("count") else "" for s in stats_list],
        [f"{s.get('mean', 0.0):.3f}" if s.get("count") else "" for s in stats_list],
        [f"{s.get('min', 0.0):.3f}" if s.get("count") else "" for s in stats_list],
        [f"{s.get('q25', 0.0):.3f}" if s.get("count") else "" for s in stats_list],
        [f"{s.get('median', 0.0):.3f}" if s.get("count") else "" for s in stats_list],
        [f"{s.get('q75', 0.0):.3f}" if s.get("count") else "" for s in stats_list],
        [f"{s.get('max', 0.0):.3f}" if s.get("count") else "" for s in stats_list],
        [
            f"{val:.3f}" if val is not None else ""
            for val in (_bowley(s) for s in stats_list)
        ],
    ]
    table = ax.table(
        cellText=cell_text,
        rowLabels=stats_rows,
        colLabels=labels,
        cellLoc="center",
        rowLoc="center",
        loc="bottom",
        bbox=[0.0, -0.6, 1.0, 0.5],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(sizes["table"])
    col_colors = _component_color_map()
    for (row_idx, col_idx), cell in table.get_celld().items():
        if col_idx in col_colors:
            cell.get_text().set_color(col_colors[col_idx])
            _apply_text_outline(cell.get_text(), linewidth=0.3)
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


FILMAN_WALLS_COMPONENTS = [
    ("filament_manifolds_only", "Filaments excl. walls", "#f2e661"),
    ("shared_filament_manifolds_walls", "Filaments and walls", "#fca50a"),
    ("walls_not_filament_manifolds", "Walls excl. filaments", "#c73e4c"),
    ("unassigned_walls_filament_manifolds", "Not filaments or walls", "#5d126e"),
]


def main() -> None:
    args = parse_args()
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
        run_polars(paths, args)
        return
    if args.engine == "stream":
        run_stream(paths, args)
        return

    scalars, values = read_inputs(paths)
    if not scalars:
        raise SystemExit("[error] No scalar columns found in topology_points.csv inputs.")
    hist_scalars = _select_hist_scalars(scalars, args)

    cat_counts: Dict[str, int] = {}
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

    with open(stats_path, "w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        box_stats_map: Dict[str, Dict[str, Dict[str, float]]] = {}
        for cat in sorted(values.keys()):
            row: Dict[str, object] = {"category": cat}
            for name in scalars:
                stats = aggregate_values(values[cat].get(name, []))
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
    print(f"[done] wrote {stats_path}")

    # Histograms
    try:
        import numpy as np
    except Exception as exc:
        raise SystemExit(f"[error] numpy is required for histograms: {exc}")

    have_plots = False
    if not args.no_plots:
        try:
            import matplotlib.pyplot as plt

            have_plots = True
        except Exception:
            have_plots = False
    write_filman_walls_composition_plot(
        out_dir,
        args.output_prefix,
        cat_counts,
        have_plots,
        args.plot_fontscale,
        args.plot_dpi,
    )

    edges_map: Dict[str, Dict[str, List[float]]] = {}
    plow = phigh = None
    if args.hist_percentile_range:
        plow, phigh = args.hist_percentile_range
        plow = plow / 100.0
        phigh = phigh / 100.0
    qlo = 0.0 if plow is None else plow
    qhi = 1.0 if phigh is None else phigh

    global_support: Dict[str, Tuple[float, float]] = {}
    for name in hist_scalars:
        all_vals: List[float] = []
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

    hist_counts: Dict[str, Dict[str, List[int]]] = {}
    hist_edges: Dict[str, Dict[str, List[float]]] = {}
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
                centers = [(edges[i] + edges[i + 1]) / 2 for i in range(len(edges) - 1)]
                widths = [edges[i + 1] - edges[i] for i in range(len(edges) - 1)]
                fig = plt.figure()
                plt.bar(centers, counts, width=widths, align="center")
                plt.axvline(float(np.median(arr)), color="red", linestyle="--", linewidth=1)
                plt.xlabel(name)
                plt.ylabel("count")
                plt.title(f"{cat}: {name} (n={len(vals):.2e})")
                fig_path = out_dir / f"{args.output_prefix}_{cat}_{name}_hist.png"
                fig.savefig(fig_path, dpi=args.plot_dpi, bbox_inches="tight")
                plt.close(fig)
    if have_plots:
        plot_scalars = {args.violin_scalar, args.box_scalar}
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
        if args.violin_scalar in plot_values:
            write_filman_walls_violin_plot(
                out_dir,
                args.output_prefix,
                plot_values[args.violin_scalar],
                have_plots,
                scalar_name=args.violin_scalar,
                font_scale=args.plot_fontscale,
                plot_dpi=args.plot_dpi,
            )
        if args.box_scalar in plot_values:
            write_filman_walls_box_plot(
                out_dir,
                args.output_prefix,
                box_stats_map.get(args.box_scalar, {}),
                have_plots,
                scalar_name=args.box_scalar,
                font_scale=args.plot_fontscale,
                plot_dpi=args.plot_dpi,
            )



if __name__ == "__main__":
    main()
