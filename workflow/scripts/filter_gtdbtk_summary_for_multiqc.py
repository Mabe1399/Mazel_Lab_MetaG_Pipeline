#!/usr/bin/env python3
# Author: Matias Becker Burgos (Matias.BeckerBurgos@unil.ch)
#
# Usage: python3 filter_gtdbtk_summary_for_multiqc.py --input gtdbtk.bac120.summary.tsv --output gtdbtk_multiqc.tsv
#
# gtdbtk.bac120.summary.tsv has ~20 columns (ANI/AF against several reference
# genomes, MSA percent, translation table, ...) — too wide for a skimmable
# MultiQC table. This keeps only the columns that actually matter for judging
# classification confidence at a glance: the call itself, how it was made,
# and anything flagged as unusual.

import argparse
import pandas as pd

KEEP_COLUMNS = ["user_genome", "classification", "classification_method", "note", "warnings"]

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

summary = pd.read_csv(args.input, sep="\t", na_values="N/A")
filtered = summary[[c for c in KEEP_COLUMNS if c in summary.columns]]
filtered.to_csv(args.output, sep="\t", index=False)
