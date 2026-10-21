# FINAL PIPELINE-WIDE MULTIQC REPORT SNAKEFILE
# Last Updated: 2026-07-28
# Author: Matias Becker Burgos (Matias.BeckerBurgos@unil.ch)

#######################################################
# DESCRIPTION
#######################################################
# Aggregates all pipeline outputs into a single MultiQC report.
#
# MODULE INTERFACE
#   INPUT:       Outputs from all enabled pipeline stages
#   OUTPUT:      {OUTPUT_DIR}/02_Results/MultiQC_FinalReport/multiqc_report.html
#   REQUIRES:    At minimum: QC: true
#   STEPS KEY:   final_report: true
#   USES:        SAMPLES, CONTIGS
#   NOTE:        Benchmark aggregation runs first inside the shell block so the
#                benchmark custom section is included in the report. Snakemake
#                input declarations are dependency anchors, not the exhaustive
#                file list — MultiQC auto-discovers reports from the directories.


import os as _os


def get_multiqc_inputs(wildcards):
    """Return input files for dependency tracking only (not the exhaustive MultiQC file list)."""
    inputs = {}

    # FastQC zips as dependency anchors for QC stage completion
    if config["steps"].get("QC", True):
        inputs["fastqc_after"] = expand(
            f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_After/{{sample}}_R{{read}}_paired_fastqc.zip",
            sample=SAMPLES,
            read=[1, 2],
        )
        # Per-sample Kraken2 reports are temp()'d once parse_kraken2_report consumes
        # them; anchoring here too defers that cleanup until this rule has also run,
        # since MultiQC's "kraken" module parses them directly from 01_Analysis/.
        inputs["kraken2_reports"] = expand(
            f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/Kraken2/{{sample}}_kraken2.report",
            sample=SAMPLES,
        )
    else:
        inputs["fastqc_check"] = expand(
            f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_Check/{{sample}}_R{{read}}_cleaned_fastqc.zip",
            sample=SAMPLES,
            read=[1, 2],
        )

    # MetaQuast runs for all assembly samples (no steps gate on assembly)
    inputs["metaquast"] = expand(
        f"{OUTPUT_DIR}/01_Analysis/02_Assembly/Metaquast_QC/{{sample}}",
        sample=CONTIGS,
    )

    # CheckM is always part of the binning pipeline
    inputs["checkm"] = f"{OUTPUT_DIR}/02_Results/03_Binning/CheckM_QC/all_MAGs_stats.tsv"
    # Per-sample qa_summary.tsv are temp()'d once Bin_filtering_QC_summary consumes
    # them; anchoring here too defers that cleanup until this rule has also run,
    # since MultiQC's "checkm" module parses them directly from 01_Analysis/.
    inputs["checkm_raw"] = expand(
        f"{OUTPUT_DIR}/01_Analysis/03_Binning/CheckM_QC/{{contig}}_checkm_QC/qa_summary.tsv",
        contig=CONTIGS,
    )

    if config["steps"].get("MAGs_Reporting_step"):
        inputs["mag_info"] = (
            f"{OUTPUT_DIR}/02_Results/04_MAGs_Reporting_Annotating/MAG_Representative_info_parsed.tsv"
        )

    if config["steps"].get("taxonomic_profiling"):
        inputs["kraken2_tax"] = (
            f"{OUTPUT_DIR}/02_Results/02b_Taxonomic_Profiling/Kraken2/"
            f"Combined_kraken2_tax_profiling_report.tsv"
        )

    return inputs


def _build_multiqc_dirs():
    """Build the list of directories/files to pass to MultiQC."""
    dirs = []
    if config["steps"].get("QC", True):
        dirs.append(f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/")
    else:
        dirs.append(f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_Check/")
    dirs += [
        f"{OUTPUT_DIR}/01_Analysis/02_Assembly/Metaquast_QC/",
        f"{OUTPUT_DIR}/01_Analysis/03_Binning/CheckM_QC/",
    ]
    if config["steps"].get("QC", True):
        dirs += [
            f"{LOG_DIR}/01_QC_Preprocessing/Trimming/",
            f"{LOG_DIR}/01_QC_Preprocessing/HF/",
        ]
    if config["steps"].get("MAGs_Reporting_step"):
        dirs.append(f"{OUTPUT_DIR}/01_Analysis/04_MAGs_Reporting_Annotating/")

    return dirs


MULTIQC_SEARCH_DIRS = _build_multiqc_dirs()


rule final_multiqc:
    input:
        unpack(get_multiqc_inputs)
    output:
        report = f"{OUTPUT_DIR}/02_Results/MultiQC_FinalReport/multiqc_report.html"
    params:
        outdir      = f"{OUTPUT_DIR}/02_Results/MultiQC_FinalReport",
        bench_dir   = BENCH_DIR,
        bench_tsv   = f"{OUTPUT_DIR}/02_Results/benchmark_summary.tsv",
        search_dirs = lambda wc: " ".join(MULTIQC_SEARCH_DIRS),
        config_file = f"{workflow.basedir}/../config/multiqc_config.yaml",
    conda:
        "../envs/QC_env.yaml"
    log:
        f"{LOG_DIR}/06_FinalReport/final_multiqc.log"
    benchmark:
        f"{BENCH_DIR}/06_FinalReport/final_multiqc.tsv"
    threads: 1
    localrule: True
    shell:
        """
        set -euo pipefail

        python {workflow.basedir}/scripts/aggregate_benchmarks.py \
            --bench-dir {params.bench_dir} \
            --output {params.bench_tsv} > {log} 2>&1

        mkdir -p {params.outdir}
        multiqc \
            {params.search_dirs} \
            {params.bench_tsv} \
            --config {params.config_file} \
            --outdir {params.outdir} \
            --filename multiqc_report.html \
            --force \
            --quiet >> {log} 2>&1
        """
