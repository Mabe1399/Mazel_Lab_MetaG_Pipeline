# CLEAN DATA CHECK SNAKEFILE
# Last Updated: 2026-06-12
# Author: Matias Becker Burgos (Matias.BeckerBurgos@unil.ch)

#######################################################
# DESCRIPTION
#######################################################
# MODULE INTERFACE
#   INPUT:    Pre-cleaned FASTQ files in CLEAN_DATA_DIR (= DATA_DIR when QC: false)
#             Expected pattern: {sample}_R{1,2}_cleaned.fastq.gz
#   OUTPUT:   {OUTPUT_DIR}/02_Results/01_QC_Preprocessing/QC_Check/multiqc_report.html
#   REQUIRES: QC: false  (mutually exclusive with 01_quality_control_preprocessing.smk)
#   STEPS KEY: check_clean_data
#   USES:     SAMPLES

# Lightweight sanity check: runs FastQC + MultiQC on pre-cleaned reads to
# confirm adapter removal, trimming quality, and absence of host contamination
# before committing to assembly.

##################################################
# RULES
##################################################

rule check_fastqc:
    input:
        lambda wildcards: [
            f"{CLEAN_DATA_DIR}/{wildcards.sample}_R1_cleaned.fastq.gz",
            f"{CLEAN_DATA_DIR}/{wildcards.sample}_R2_cleaned.fastq.gz",
        ]
    output:
        html1 = temp(f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_Check/{{sample}}_R1_cleaned_fastqc.html"),
        html2 = temp(f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_Check/{{sample}}_R2_cleaned_fastqc.html"),
        zip1  = f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_Check/{{sample}}_R1_cleaned_fastqc.zip",
        zip2  = f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_Check/{{sample}}_R2_cleaned_fastqc.zip",
    conda:
        "../envs/QC_env.yaml"
    threads: config["fastqc"]["threads"]
    log:
        f"{LOG_DIR}/01_QC_Preprocessing/QC_Check/{{sample}}_fastqc.log"
    benchmark:
        f"{BENCH_DIR}/01_QC_Preprocessing/QC_Check/{{sample}}_fastqc.tsv"
    resources:
        mem_mb        = config["fastqc"]["memory_mb"],
        runtime       = config["fastqc"]["runtime_min"],
        cpus_per_task = config["fastqc"]["threads"]
    shell:
        """
        mkdir -p {OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_Check
        fastqc -t {threads} \
            -o {OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_Check \
            {input} > {log} 2>&1
        """


rule check_multiqc:
    input:
        expand(
            f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_Check/{{sample}}_R{{read}}_cleaned_fastqc.html",
            sample=SAMPLES,
            read=[1, 2],
        )
    output:
        outdir = directory(f"{OUTPUT_DIR}/02_Results/01_QC_Preprocessing/QC_Check")
    conda:
        "../envs/QC_env.yaml"
    threads: config["multiqc"]["threads"]
    localrule: True
    log:
        f"{LOG_DIR}/01_QC_Preprocessing/QC_Check/multiqc.log"
    benchmark:
        f"{BENCH_DIR}/01_QC_Preprocessing/QC_Check/multiqc.tsv"
    resources:
        mem_mb        = config["multiqc"]["memory_mb"],
        runtime       = config["multiqc"]["runtime_min"],
        cpus_per_task = config["multiqc"]["threads"]
    shell:
        """
        mkdir -p {output.outdir}
        multiqc \
            -o {output.outdir} \
            {OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_Check > {log} 2>&1
        """


# Mirrors 01_quality_control_preprocessing.smk's rule of the same name — the
# two files are mutually exclusive via the QC steps flag, so no collision risk,
# but _get_tool_version_inputs() in the Snakefile requires qc_env_versions.tsv
# unconditionally regardless of which branch is active.
rule tool_version_qc_env:
    output:
        versions = temp(f"{OUTPUT_DIR}/01_Analysis/ToolVersions/qc_env_versions.tsv")
    conda:
        "../envs/QC_env.yaml"
    localrule: True
    shell:
        """
        mkdir -p $(dirname {output.versions})
        {{
            printf "fastqc\\t%s\\n" "$(fastqc --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+(\\.[0-9]+)?' | head -1 || echo unknown)"
            printf "multiqc\\t%s\\n" "$(multiqc --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+(\\.[0-9]+)?' | head -1 || echo unknown)"
        }} > {output.versions}
        """
