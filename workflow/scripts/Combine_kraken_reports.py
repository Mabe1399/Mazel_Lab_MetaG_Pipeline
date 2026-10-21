#!/usr/bin/env python3
# Author: Matias Becker Burgos (Matias.BeckerBurgos@unil.ch)
# Script to combine multiple Kraken2 reports into a single report
# Usage: python Combine_kraken_reports.py -i report1.txt report2.txt ... -o combined_report.txt

# output format: taxid name rank sample1_pct sample1_taxa_reads sample1_clade_reads sample2_pct sample2_taxa_reads sample2_clade_reads ...

import argparse
import pandas as pd
from pathlib import Path


def parse_kraken_reports(file_path):
    sample_name = Path(file_path).stem.replace('_kraken2', '')

    df = pd.read_csv(
        file_path, 
        sep='\t', 
        header=None, 
        names=[
        'pct', 'reads_clade', 'reads_taxon', 'rank', 'taxid', 'name'],
        dtype={'taxid': str}
        )
    
    df = df[['taxid', 'name', 'rank', 'pct', 'reads_taxon', 'reads_clade']]
    df = df.rename(columns={
        "pct": f'{sample_name}_pct',
        "reads_taxon":  f'{sample_name}_taxa_reads',
        "reads_clade": f'{sample_name}_clade_reads'
        })
    
    return df

def main():
    parser = argparse.ArgumentParser(description='Combine multiple Kraken2 reports into a single report.')
    parser.add_argument('-i', '--input', nargs='+', required=True, help='Input Kraken2 report files.')
    parser.add_argument('-o', '--output', required=True, help='Output combined report file.')
    args = parser.parse_args()

    combined_df = pd.DataFrame()

    for file_path in args.input:
        df = parse_kraken_reports(file_path)

        if combined_df.empty:
            combined_df = df
        else:
            combined_df = pd.merge(combined_df, df, on=['taxid', 'name', 'rank'], how='outer')

    combined_df = combined_df.sort_values(by='taxid', key=lambda x: x.astype(int))

    combined_df.to_csv(args.output, sep='\t', index=False)

if __name__ == "__main__":
    main()