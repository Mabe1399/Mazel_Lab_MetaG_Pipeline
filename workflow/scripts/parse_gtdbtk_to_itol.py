#!/usr/bin/env python3
# Author: Matias Becker Burgos (Matias.BeckerBurgos@unil.ch)
# 
# Usage: python parse_gtdbtk_to_itol.py --summary gtdbtk_classify_output --metadata metadata_file \
#            --col-host Host_species --col-sample-id Glomus_ID --split-criteria _DASTool_ \
#            --host-colors '{"Apodemus ilex": "#ffe6a7"}' --outdir output_itol_dir
#
# Column names, split string and host colours come from the run's config (DATA_VARIABLES and
# host_colors), passed in by the parse_gtdbtk_summary_to_iTOL rule.

# Output file: leaf_genus_gtdb_to_iTOL_Phylum_strip_labels.txt
# Output file: leaf_genus_gtdb_to_iTOL_Genus_strip_colors.txt
# Output file: leaf_genus_gtdb_to_iTOL_Host_strip_labels.txt


import pandas as pd
import json
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import colorsys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--metadata", required=True)
parser.add_argument("--summary", required=True)
parser.add_argument("--col-host", required=True, help="metadata column with the host species")
parser.add_argument("--col-sample-id", required=True, help="metadata column with the sample IDs")
parser.add_argument("--split-criteria", required=True, help="string splitting user_genome into sample ID and bin")
parser.add_argument("--host-colors", default="{}", help="JSON object mapping host species to colours")
parser.add_argument("--outdir", required=True)
args = parser.parse_args()

if args.metadata.endswith(".xlsx"):
    metadata = pd.read_excel(args.metadata)
else:
    metadata = pd.read_csv(args.metadata, sep="\t")

gtdb_summary = pd.read_csv(args.summary, sep="\t")

HOST_COLUMN_NAME = args.col_host
SPLIT_CRITERIA = args.split_criteria
SAMPLE_ID_COL = args.col_sample_id

HOST_COLORS = json.loads(args.host_colors)

for col in (HOST_COLUMN_NAME, SAMPLE_ID_COL):
    if col not in metadata.columns:
        raise SystemExit(f"Column '{col}' not found in metadata {args.metadata}; "
                         f"check DATA_VARIABLES in the config. Available: {', '.join(metadata.columns)}")

#########################################################

gtdb_summary["Sample_ID"] = (
    gtdb_summary["user_genome"].str.split(SPLIT_CRITERIA).str[0]
    )

merged = pd.merge(
    gtdb_summary,
    metadata,
    right_on=SAMPLE_ID_COL,
    left_on="Sample_ID",
    how="left"
    )

def parse_tax(t):
    tax_levels = t.split(";")
    tax = {}
    for level in tax_levels:
        if "__" in level:
            key, value = level.split("__", 1)
            tax[key] = value
    return tax

tax_table = merged["classification"].apply(parse_tax).apply(pd.Series)
tax_table = tax_table.replace("", pd.NA)

## Leaf Genus (fall back to family if genus is NaN)
merged["leaf_genus"] = "g__" + tax_table["g"].astype(str)
merged["leaf_genus"] = merged["leaf_genus"].fillna(
    "f__" + tax_table["f"].fillna("Unclassified").astype(str)
    )

merged["leaf_phylum"] = "p__" + tax_table["p"].astype(str)

genome_to_genus = dict(zip(merged["user_genome"], merged["leaf_genus"]))
genome_to_phylum = dict(zip(merged["user_genome"], merged["leaf_phylum"]))
genome_to_host = dict(zip(merged["user_genome"], merged[HOST_COLUMN_NAME]))

genus_counts = merged["leaf_genus"].value_counts()
phylum_counts = merged["leaf_phylum"].value_counts()
host_counts = merged[HOST_COLUMN_NAME].value_counts()

leaf_id_to_Genus_label = {}
for genome, genus in genome_to_genus.items():
    if genus.startswith("g__"):
        name = genus[3:]
    elif genus.startswith("f__"):
        name = f"Unclassified_{genus[3:]}"
    else:
        name = genus

    count = genus_counts.get(genus, 0)
    leaf_id_to_Genus_label[genome] = f"{name} (n={count})"


leaf_id_to_Phylum_label = {}
for genome, phylum in genome_to_phylum.items():
    name = phylum[3:] if phylum.startswith("p__") else phylum
    count = phylum_counts.get(phylum, 0)
    leaf_id_to_Phylum_label[genome] = f"{name} (n={count})"


leaf_id_to_Host_label = {}
for genome, host in genome_to_host.items():
    host_clean = str(host).replace("_", " ")
    count = host_counts.get(host, 0)
    leaf_id_to_Host_label[genome] = f"{host_clean} (n={count})"

# --- Phylum base colors ---
phyla = sorted(merged["leaf_phylum"].unique())
tab20 = plt.cm.tab20.colors

phylum_to_color = {}
for i, phylum in enumerate(phyla):
    base_color = tab20[i % len(tab20)]
    phylum_to_color[phylum] = mcolors.to_hex(base_color)

# --- Genus shades per phylum ---
leaf_id_to_genus_and_color = {}
leaf_id_to_phylum_and_color = {}
leaf_id_to_host_and_color = {}

for phylum in phyla:
    base_rgb = mcolors.to_rgb(phylum_to_color[phylum])
    h, l, s = colorsys.rgb_to_hls(*base_rgb)

    genera_in_phylum = sorted(
        merged.loc[merged["leaf_phylum"] == phylum, "leaf_genus"].unique()
    )

    n = len(genera_in_phylum)

    for i, genus in enumerate(genera_in_phylum):
        # vary lightness within safe range
        lightness = 0.35 + (0.4 * i / max(n - 1, 1))
        r, g, b = colorsys.hls_to_rgb(h, lightness, s)
        genus_hex = mcolors.to_hex((r, g, b))

        genomes = merged.loc[
            (merged["leaf_phylum"] == phylum) &
            (merged["leaf_genus"] == genus),
            "user_genome"
        ]

        for genome in genomes:
            leaf_id_to_genus_and_color[genome] = (
                leaf_id_to_Genus_label[genome],
                genus_hex
            )

# --- Phylum strip ---
for genome in merged["user_genome"]:
    phylum = genome_to_phylum[genome]
    leaf_id_to_phylum_and_color[genome] = (
        leaf_id_to_Phylum_label[genome],
        phylum_to_color[phylum]
    )

# --- Host strip ---
for genome in merged["user_genome"]:
    host = genome_to_host[genome]
    color = HOST_COLORS.get(host, "#D3D3D3")
    leaf_id_to_host_and_color[genome] = (
        leaf_id_to_Host_label[genome],
        color
    )

###########################################################
# CREATE OUTPUT FILES
Path(args.outdir).mkdir(parents=True, exist_ok=True)

def write_itol_strip(outfile, label, branch_color, data_dict, add_legend=True):
    """
    Writes iTOL strip with optional embedded legend.
    Removes (n=...) from legend labels only.
    """

    with open(outfile, "w") as f:
        f.write("DATASET_COLORSTRIP\n")
        f.write("SEPARATOR TAB\n")
        f.write(f"DATASET_LABEL\t{label}\n")
        f.write("COLOR\t#ff0000\n")
        f.write(f"COLOR_BRANCHES\t{branch_color}\n\n")

        # ---- OPTIONAL LEGEND BLOCK ----
        if add_legend:
            legend = {}

            for _, (text, color) in data_dict.items():
                clean_label = text.split(" (n=")[0]
                legend[clean_label] = color

            legend_labels = sorted(legend.keys())
            legend_colors = [legend[l] for l in legend_labels]
            legend_shapes = ["1"] * len(legend_labels)

            f.write(f"LEGEND_TITLE\t{label}\n")
            f.write("LEGEND_SHAPES\t" + "\t".join(legend_shapes) + "\n")
            f.write("LEGEND_COLORS\t" + "\t".join(legend_colors) + "\n")
            f.write("LEGEND_LABELS\t" + "\t".join(legend_labels) + "\n\n")

        # ---- DATA BLOCK ----
        f.write("DATA\n")
        for leaf_id, (text, color) in data_dict.items():
            f.write(f"{leaf_id}\t{color}\t{text}\n")


write_itol_strip(
    Path(args.outdir) / "leaf_phylum_gtdb_to_iTOL_strip_labels.txt",
    "Leaf_Phylum",
    1,
    leaf_id_to_phylum_and_color,
    add_legend=True
)

write_itol_strip(
    Path(args.outdir) / "leaf_genus_gtdb_to_iTOL_strip_labels.txt",
    "Leaf_Genus",
    1,
    leaf_id_to_genus_and_color,
    add_legend=False
)

write_itol_strip(
    Path(args.outdir) / "leaf_host_gtdb_to_iTOL_strip_labels.txt",
    "Host_Species",
    0,
    leaf_id_to_host_and_color,
    add_legend=True
)