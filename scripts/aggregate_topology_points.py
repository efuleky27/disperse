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
from typing import Dict, Iterable, List, Tuple


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Aggregate topology_points.csv files across crops.")
    p.add_argument("--root", default="outputs", help="Root directory to search for topology_points.csv files.")
    p.add_argument("--glob", default="**/*_topology_points.csv", help="Glob pattern under --root.")
    p.add_argument("--inputs", nargs="+", help="Explicit list of topology_points.csv files.")
    p.add_argument("--output-dir", required=True, help="Directory for combined outputs.")
    p.add_argument("--output-prefix", default="topology_combined", help="Prefix for output files.")
    p.add_argument("--bins", type=int, default=100, help="Histogram bins per scalar.")
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
    p.add_argument("--no-plots", action="store_true", help="Skip PNG histogram plots.")
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


def read_inputs(paths: Iterable[Path]) -> Tuple[List[str], Dict[str, Dict[str, List[float]]]]:
    scalars: List[str] = []
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

    if args.no_plots:
        have_plots = False
    else:
        try:
            import matplotlib.pyplot as plt  # noqa: F401

            have_plots = True
        except Exception:
            have_plots = False

    for name in scalars:
        if args.hist_percentile_range:
            plow, phigh = args.hist_percentile_range
            bounds = data.group_by("category").agg(
                [
                    pl.col(name).quantile(plow / 100.0, "linear").alias("low"),
                    pl.col(name).quantile(phigh / 100.0, "linear").alias("high"),
                ]
            )
        else:
            bounds = data.group_by("category").agg(
                [pl.col(name).min().alias("low"), pl.col(name).max().alias("high")]
            )
        bounds_df = _collect(bounds)
        bounds_map = {row["category"]: (row["low"], row["high"]) for row in bounds_df.to_dicts()}
        bounds_lazy = bounds
        width = pl.col("high") - pl.col("low")
        bin_expr = (
            ((pl.col(name) - pl.col("low")) / width * args.bins)
            .floor()
            .cast(pl.Int64)
            .clip(0, args.bins - 1)
            .alias("bin")
        )
        hist_df = (
            data.join(bounds_lazy, on="category", how="left")
            .filter(width > 0)
            .with_columns(bin_expr)
            .group_by(["category", "bin"])
            .agg(pl.count().alias("count"))
        )
        hist_df = _collect(hist_df)
        for cat, (lo, hi) in bounds_map.items():
            if hi is None or lo is None or hi <= lo:
                continue
            counts = [0] * args.bins
            subset = hist_df.filter(pl.col("category") == cat)
            for row in subset.to_dicts():
                counts[int(row["bin"])] = int(row["count"])
            edges = [lo + (hi - lo) * i / args.bins for i in range(args.bins + 1)]
            csv_path = out_dir / f"{args.output_prefix}_{cat}_{name}_hist.csv"
            write_histogram_csv(csv_path, edges, counts)
            if have_plots:
                import matplotlib.pyplot as plt

                fig = plt.figure()
                centers = [(edges[i] + edges[i + 1]) / 2 for i in range(len(edges) - 1)]
                widths = [edges[i + 1] - edges[i] for i in range(len(edges) - 1)]
                plt.bar(centers, counts, width=widths, align="center")
                plt.xlabel(name)
                plt.ylabel("count")
                plt.title(f"{cat}: {name}")
                fig_path = out_dir / f"{args.output_prefix}_{cat}_{name}_hist.png"
                fig.savefig(fig_path, dpi=150, bbox_inches="tight")
                plt.close(fig)


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
        if is_wall and is_filman:
            cats.append("shared_walls_filament_manifolds")
        if is_fil and is_filman:
            cats.append("shared_filaments_filament_manifolds")
        if is_wall and not is_filman:
            cats.append("walls_not_filament_manifolds")
        if is_filman and not is_wall:
            cats.append("filament_manifolds_not_walls")
    if not is_wall and not is_fil and (not include_filman or not is_filman):
        cats.append("unassigned")
    return cats


def run_stream(paths: List[Path], args: argparse.Namespace) -> None:
    import math
    import random

    sample_size = 20000
    scalars: List[str] = []
    include_filman = False
    stats: Dict[str, Dict[str, Dict[str, float]]] = {}
    category_count: Dict[str, int] = {}
    samples: Dict[str, Dict[str, List[float]]] = {}

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
                        if args.hist_percentile_range:
                            samp = _sample_store(cat, name)
                            seen = int(store["count"])
                            if len(samp) < sample_size:
                                samp.append(val)
                            else:
                                j = random.randint(1, seen)
                                if j <= sample_size:
                                    samp[j - 1] = val

    if not scalars:
        raise SystemExit("[error] No scalar columns found in topology_points.csv inputs.")

    # histogram bounds
    bounds: Dict[str, Dict[str, Tuple[float, float]]] = {}
    for cat, scal_map in stats.items():
        for name, store in scal_map.items():
            if store["count"] <= 0:
                continue
            low = store["min"]
            high = store["max"]
            if args.hist_percentile_range:
                samp = samples.get(cat, {}).get(name, [])
                if samp:
                    samp_sorted = sorted(samp)
                    n = len(samp_sorted)
                    plow, phigh = args.hist_percentile_range
                    lo_idx = int((plow / 100.0) * (n - 1))
                    hi_idx = int((phigh / 100.0) * (n - 1))
                    low = samp_sorted[lo_idx]
                    high = samp_sorted[hi_idx]
            if high > low:
                bounds.setdefault(cat, {})[name] = (low, high)

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
                    for cat in cats:
                        if name not in bounds.get(cat, {}):
                            continue
                        low, high = bounds[cat][name]
                        if val < low or val > high:
                            continue
                        width = high - low
                        if width <= 0:
                            continue
                        bin_idx = int((val - low) / width * args.bins)
                        if bin_idx >= args.bins:
                            bin_idx = args.bins - 1
                        counts.setdefault(cat, {}).setdefault(name, [0] * args.bins)
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
                # approximate quantiles from histogram bins
                qvals = {"q25": 0.0, "median": 0.0, "q75": 0.0}
                if name in counts.get(cat, {}) and name in bounds.get(cat, {}):
                    hist = counts[cat][name]
                    total = sum(hist)
                    if total > 0:
                        low, high = bounds[cat][name]
                        edges = [low + (high - low) * i / args.bins for i in range(args.bins + 1)]
                        targets = {
                            "q25": 0.25 * total,
                            "median": 0.5 * total,
                            "q75": 0.75 * total,
                        }
                        cum = 0
                        for i, c in enumerate(hist):
                            prev = cum
                            cum += c
                            for key, target in list(targets.items()):
                                if prev <= target <= cum and c > 0:
                                    frac = (target - prev) / c
                                    qvals[key] = edges[i] + frac * (edges[i + 1] - edges[i])
                                    targets.pop(key, None)
                            if not targets:
                                break
                row[f"{name}_q25"] = qvals["q25"]
                row[f"{name}_median"] = qvals["median"]
                row[f"{name}_q75"] = qvals["q75"]
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

    for cat, name_map in counts.items():
        for name, hist in name_map.items():
            if name not in bounds.get(cat, {}):
                continue
            low, high = bounds[cat][name]
            edges = [low + (high - low) * i / args.bins for i in range(args.bins + 1)]
            csv_path = out_dir / f"{args.output_prefix}_{cat}_{name}_hist.csv"
            write_histogram_csv(csv_path, edges, hist)
            if have_plots:
                import matplotlib.pyplot as plt

                centers = [(edges[i] + edges[i + 1]) / 2 for i in range(len(edges) - 1)]
                widths = [edges[i + 1] - edges[i] for i in range(len(edges) - 1)]
                fig = plt.figure()
                plt.bar(centers, hist, width=widths, align="center")
                plt.xlabel(name)
                plt.ylabel("count")
                plt.title(f"{cat}: {name}")
                fig_path = out_dir / f"{args.output_prefix}_{cat}_{name}_hist.png"
                fig.savefig(fig_path, dpi=150, bbox_inches="tight")
                plt.close(fig)


def write_histogram_csv(path: Path, edges: List[float], counts: List[int]) -> None:
    with open(path, "w", newline="", encoding="ascii") as handle:
        writer = csv.writer(handle)
        writer.writerow(["bin_left", "bin_right", "count"])
        for i, count in enumerate(counts):
            writer.writerow([edges[i], edges[i + 1], count])


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

    for cat, cat_vals in values.items():
        for name in scalars:
            vals = cat_vals.get(name, [])
            if not vals:
                continue
            arr = np.asarray(vals, dtype=float)
            if args.hist_percentile_range:
                lo, hi = np.percentile(arr, args.hist_percentile_range)
            else:
                lo, hi = float(arr.min()), float(arr.max())
            if hi <= lo:
                continue
            counts, edges = np.histogram(arr, bins=args.bins, range=(lo, hi))
            csv_path = out_dir / f"{args.output_prefix}_{cat}_{name}_hist.csv"
            write_histogram_csv(csv_path, edges.tolist(), counts.tolist())
            if have_plots:
                fig = plt.figure()
                plt.hist(arr, bins=args.bins, range=(lo, hi))
                plt.xlabel(name)
                plt.ylabel("count")
                plt.title(f"{cat}: {name}")
                fig_path = out_dir / f"{args.output_prefix}_{cat}_{name}_hist.png"
                fig.savefig(fig_path, dpi=150, bbox_inches="tight")
                plt.close(fig)


if __name__ == "__main__":
    main()
