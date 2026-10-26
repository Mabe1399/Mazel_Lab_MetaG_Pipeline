#!/usr/bin/env python3
# Author: Matias Becker Burgos (Matias.BeckerBurgos@unil.ch)
#
# Usage: python3 aggregate_dastool_summary.py --input sample1_DASTool_summary.tsv sample2_DASTool_summary.tsv ... \
#            --samples sample1 sample2 ... --output DASTool_binning_summary.tsv
#
# Each {contig}_DASTool_summary.tsv already lists only the bins DAS_Tool actually
# selected into the consensus (one row per final bin, not every candidate from
# every binner — confirmed against DAS_Tool 1.1.7's own source), with a bin_set
# column naming which binner(s) it came from. This just concatenates those
# per-sample tables with a sample column and keeps the columns worth showing
# in a quick-skim report.

import argparse
import pandas as pd

KEEP_COLUMNS = ["bin", "bin_set", "bin_score", "SCG_completeness", "SCG_redundancy"]

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True, nargs="+")
parser.add_argument("--samples", required=True, nargs="+")
parser.add_argument("--output", required=True)
args = parser.parse_args()

frames = []
for sample, path in zip(args.samples, args.input):
    df = pd.read_csv(path, sep="\t")
    df = df[[c for c in KEEP_COLUMNS if c in df.columns]]
    df.insert(0, "sample", sample)
    frames.append(df)

combined = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(
    columns=["sample"] + KEEP_COLUMNS
)
combined.to_csv(args.output, sep="\t", index=False)
