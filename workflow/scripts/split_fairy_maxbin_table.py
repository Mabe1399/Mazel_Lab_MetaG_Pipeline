#!/usr/bin/env python3
# Author: Matias Becker Burgos (Matias.BeckerBurgos@unil.ch)
#
# Usage: python3 split_fairy_maxbin_table.py --input wide_table.tsv --outdir coverage_dir \
#            --samples sample1 sample2 ... --output-list abund_list.txt
#
# fairy coverage --maxbin-format produces a single wide table (contigName, sample1,
# sample2, ...); MaxBin2's -abund_list needs one 2-column (contig, coverage) file per
# sample instead. This splits the wide table by column name (not position, since
# column order isn't guaranteed to match --samples) and writes the list of resulting
# file paths MaxBin2 expects.

import argparse
import os
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--outdir", required=True)
parser.add_argument("--samples", required=True, nargs="+")
parser.add_argument("--output-list", required=True)
args = parser.parse_args()

wide = pd.read_csv(args.input, sep="\t")

os.makedirs(args.outdir, exist_ok=True)

with open(args.output_list, "w") as list_f:
    for sample in args.samples:
        per_sample = wide[["contigName", sample]]
        out_path = os.path.join(args.outdir, f"{sample}_coverage.txt")
        per_sample.to_csv(out_path, sep="\t", index=False, header=False)
        list_f.write(f"{out_path}\n")
