#!/usr/bin/env python3
"""
aggregate_benchmarks.py — Summarize Snakemake benchmark files across a run.

Walks BENCH_DIR, reads all .tsv benchmark files, and produces a consolidated
summary table grouped by pipeline module. Flags the top-5 most expensive rules
by wall time and peak memory.

Usage:
    python aggregate_benchmarks.py --bench-dir /path/to/benchmark --output summary.tsv
"""

import argparse
import os
import sys
import pandas as pd
from pathlib import Path


BENCHMARK_COLS = ["s", "h:m:s", "max_rss", "max_vms", "max_uss", "max_pss",
                  "io_in", "io_out", "mean_load", "cpu_time"]


def parse_args():
    p = argparse.ArgumentParser(description="Aggregate Snakemake benchmark files.")
    p.add_argument("--bench-dir", required=True, help="Root benchmark directory")
    p.add_argument("--output",    required=True, help="Output TSV path")
    p.add_argument("--top-n",     type=int, default=5, help="Number of top bottlenecks to flag")
    return p.parse_args()


def read_benchmark(path: Path) -> dict | None:
    try:
        df = pd.read_csv(path, sep="\t")
        if df.empty or "s" not in df.columns:
            return None
        row = df.iloc[0]
        return {
            "wall_time_s":    float(row.get("s", 0)),
            "max_rss_mb":     float(row.get("max_rss", 0)),
            "cpu_time_s":     float(row.get("cpu_time", 0)),
        }
    except Exception:
        return None


def infer_module_rule_sample(rel_path: str) -> tuple[str, str, str]:
    """Extract (module, rule, sample) from relative benchmark path."""
    parts = Path(rel_path).parts
    # Expected structure: module/rule/sample.tsv  or  module/rule.tsv
    module = parts[0] if len(parts) >= 1 else "unknown"
    rule   = parts[1] if len(parts) >= 2 else "unknown"
    sample = Path(parts[2]).stem if len(parts) >= 3 else Path(parts[-1]).stem
    return module, rule, sample


def main():
    args = parse_args()
    bench_root = Path(args.bench_dir)

    if not bench_root.exists():
        sys.exit(f"Benchmark directory not found: {bench_root}")

    records = []
    for tsv in bench_root.rglob("*.tsv"):
        rel = tsv.relative_to(bench_root)
        data = read_benchmark(tsv)
        if data is None:
            continue
        module, rule, sample = infer_module_rule_sample(str(rel))
        records.append({
            "module":         module,
            "rule":           rule,
            "sample":         sample,
            "wall_time_min":  round(data["wall_time_s"] / 60, 2),
            "max_memory_gb":  round(data["max_rss_mb"] / 1024, 3),
            "cpu_time_min":   round(data["cpu_time_s"] / 60, 2),
        })

    if not records:
        print("No benchmark files found — nothing to summarize.", file=sys.stderr)
        pd.DataFrame(columns=["module","rule","sample","wall_time_min","max_memory_gb","cpu_time_min"]).to_csv(
            args.output, sep="\t", index=False
        )
        return

    df = pd.DataFrame(records).sort_values(["module", "rule", "sample"])
    df.to_csv(args.output, sep="\t", index=False)

    module_totals = (
        df.groupby("module")[["wall_time_min", "max_memory_gb", "cpu_time_min"]]
        .agg({"wall_time_min": "sum", "max_memory_gb": "max", "cpu_time_min": "sum"})
        .rename(columns={"wall_time_min": "total_wall_min", "max_memory_gb": "peak_memory_gb", "cpu_time_min": "total_cpu_min"})
        .sort_values("total_wall_min", ascending=False)
    )

    top_wall = (
        df.groupby(["module", "rule"])["wall_time_min"]
        .max()
        .sort_values(ascending=False)
        .head(args.top_n)
    )
    top_mem = (
        df.groupby(["module", "rule"])["max_memory_gb"]
        .max()
        .sort_values(ascending=False)
        .head(args.top_n)
    )

    print(f"\nBenchmark summary written to: {args.output}")
    print(f"Total rules benchmarked: {len(df)}")
    print(f"\n{'='*60}")
    print("PER-MODULE TOTALS")
    print(f"{'='*60}")
    print(module_totals.to_string())
    print(f"\n{'='*60}")
    print(f"TOP {args.top_n} SLOWEST RULES (max wall time)")
    print(f"{'='*60}")
    for (mod, rule), val in top_wall.items():
        print(f"  {mod}/{rule}: {val:.1f} min")
    print(f"\n{'='*60}")
    print(f"TOP {args.top_n} MOST MEMORY-HUNGRY RULES (peak GB)")
    print(f"{'='*60}")
    for (mod, rule), val in top_mem.items():
        print(f"  {mod}/{rule}: {val:.2f} GB")
    print()


if __name__ == "__main__":
    main()
