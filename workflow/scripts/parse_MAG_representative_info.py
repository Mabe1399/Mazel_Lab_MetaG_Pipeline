#!/usr/bin/env python3
# Author: Matias Becker Burgos (Matias.BeckerBurgos@unil.ch)
# 
# Usage: python3 parse_MAG_representative_info.py --CheckM_summary CheckM_summary_file.tsv \
# --GTDB_Tk_summary GTDB_Tk_summary_file.tsv \
# --drep_dir drep_directory \
# --output output_file.tsv

# Output file: MAG_representative_info_parsed.tsv
# This script parses the CheckM2 summary, GTDB-Tk summary and dRep directory to extract information about the representative MAGs and others, including their genome ID, completeness, contamination, taxonomy and dRep cluster.
# The output file will contain the following columns: MAG ID, completeness, contamination, Completeness_Model_Used, Domain, Phylum, Class, Order, Family, Genus, Species, dRep_Score, dRep_cluster, closest_cluster_member, furthest_cluster_member.

import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("--CheckM_summary", required=True)
parser.add_argument("--GTDB_Tk_summary", required=True)
parser.add_argument("--drep_dir", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

checkm_summary = pd.read_csv(args.CheckM_summary, sep="\t")

gtdbtk_summary = pd.read_csv(args.GTDB_Tk_summary, sep="\t", na_values="N/A")

genome_info = checkm_summary[["Name", "Completeness", "Contamination", "Completeness_Model_Used" ]]
genome_info.columns = ["MAG ID", "completeness", "contamination", "Completeness_Model_Used"]
genome_info.loc[:, "MAG ID"] = genome_info["MAG ID"].apply(lambda x: x + ".fa")

taxonomy_info = gtdbtk_summary[["user_genome", "classification"]]
taxonomy_info.columns = ["MAG ID", "taxonomy"]
taxonomy_info.loc[:, "MAG ID"] = taxonomy_info["MAG ID"].apply(lambda x: x + ".fa")
taxonomy_info[["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]] = taxonomy_info["taxonomy"].str.split(";", expand=True)
taxonomy_info = taxonomy_info.drop(columns=["taxonomy"])

drep_data_tables_dir = os.path.join(args.drep_dir, "data_tables")
sdb_file = os.path.join(drep_data_tables_dir, "Sdb.csv")
sdb_df = pd.read_csv(sdb_file)
drep_cluster_info = sdb_df[["genome", "score"]]
drep_cluster_info.columns = ["MAG ID", "dRep_Score"]

cdb_file = os.path.join(drep_data_tables_dir, "Cdb.csv")
cdb_df = pd.read_csv(cdb_file)
drep_cluster_info_cdb = cdb_df[["genome", "secondary_cluster"]]
drep_cluster_info_cdb.columns = ["MAG ID", "dRep_cluster"]

widb_file = os.path.join(drep_data_tables_dir, "Widb.csv")
widb_df = pd.read_csv(widb_file)
drep_cluster_info_widb = widb_df[["genome","closest_cluster_member", "furthest_cluster_member"]]
drep_cluster_info_widb.columns = ["MAG ID","closest_cluster_member", "furthest_cluster_member"]

# MERGE ALL DATAFRAMES if info NoT available for a MAG, fill with NA
merged_df = pd.merge(genome_info, taxonomy_info, on="MAG ID", how="left")
merged_df = pd.merge(merged_df, drep_cluster_info, on="MAG ID", how="left")
merged_df = pd.merge(merged_df, drep_cluster_info_cdb, on="MAG ID", how="left")
merged_df = pd.merge(merged_df, drep_cluster_info_widb, on="MAG ID", how="left")

merged_df.to_csv(args.output, sep="\t", index=False)