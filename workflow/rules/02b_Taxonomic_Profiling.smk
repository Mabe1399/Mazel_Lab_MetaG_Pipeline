# TAXONOMIC PROFILING SNAKEFILE
# Last Updated: 2026-06-06
# Author: Matias Becker Burgos (Matias.BeckerBurgos@unil.ch)

#######################################################
# DESCRIPTION
#######################################################
# Taxonomic profiling rules usable in many different pipelines.
# Assumptions:
#   - Cleaned FASTQ files are in CLEAN_DATA_DIR (OUTPUT_DIR/00_Data/02_Clean_Data when QC runs,
#     DATA_DIR otherwise)
#   - A Kraken2 database is available at config[kraken2_db][DB_dir]
#
# MODULE INTERFACE
#   INPUT:       Cleaned FASTQ reads from 01_quality_control_preprocessing.smk
#   OUTPUT:      {OUTPUT_DIR}/02_Results/01_QC_Preprocessing/Kraken2/Combined_kraken2_report.tsv
#   REQUIRES:    QC: true (01_quality_control_preprocessing.smk)
#   STEPS KEY:   taxonomic_profiling: true
#   USES:        SAMPLES (all samples including blanks)

# Steps included:
# 01: Taxonomic profiling (Kraken2)


##################################################
# RULES
##################################################

################ 01 TAXONOMIC PROFILING ##########

rule kraken2_Taxonomic_Profiling:
    input:
        R1=f"{CLEAN_DATA_DIR}/{{sample}}_R1_cleaned.fastq.gz",
        R2=f"{CLEAN_DATA_DIR}/{{sample}}_R2_cleaned.fastq.gz",
        db= config["kraken2_db"]["DB_dir"]
    output:
        tab=temp(f"{SCRATCH_DIR}/Kraken2/{{sample}}_kraken2.out"),
        rep=temp(f"{OUTPUT_DIR}/01_Analysis/02b_Taxonomic_Profiling/Kraken2/{{sample}}_kraken2.report")
    conda:
        "../envs/Kraken2.yaml"
    threads: config["kraken2_run"]["threads"]
    log:
        f"{LOG_DIR}/02b_Taxonomic_Profiling/Kraken2/{{sample}}_kraken2_run.log"
    benchmark:
        f"{BENCH_DIR}/02b_Taxonomic_Profiling/Kraken2/{{sample}}_kraken2_run.tsv"
    resources:
        mem_mb = config["kraken2_run"]["memory_mb"],
        runtime = config["kraken2_run"]["runtime_min"],
        cpus_per_task = config["kraken2_run"]["threads"]
    shell:
        "kraken2 --use-names --threads {threads} \
        --db {input.db} \
        --report {output.rep} --gzip-compressed \
        --report-zero-counts \
        --paired {input.R1} {input.R2} \
        --output {output.tab} > {log} 2>&1"


rule parse_kraken2_taxonomic_report:
    input:
        expand(f"{OUTPUT_DIR}/01_Analysis/02b_Taxonomic_Profiling/Kraken2/{{sample}}_kraken2.report", sample=SAMPLES)
    output:
        f"{OUTPUT_DIR}/02_Results/02b_Taxonomic_Profiling/Kraken2/Combined_kraken2_tax_profiling_report.tsv"
    localrule: True
    conda:
        "../envs/Parsing.yaml"
    threads: 1
    log:
        f"{LOG_DIR}/02b_Taxonomic_Profiling/Kraken2/Combined_kraken2_parsing.log"
    benchmark:
        f"{BENCH_DIR}/02b_Taxonomic_Profiling/Kraken2/Combined_kraken2_parsing.tsv"
    shell:
        "python3 scripts/Combine_kraken_reports.py -o {output} -i {input} > {log} 2>&1"


# Defined here unconditionally (this file only loads when taxonomic_profiling:
# true) — see 01_quality_control_preprocessing.smk's matching guard for why
# the two files can't both define this rule when QC and taxonomic_profiling
# are both enabled.
rule tool_version_kraken2:
    output:
        versions = temp(f"{OUTPUT_DIR}/01_Analysis/ToolVersions/kraken2_versions.tsv")
    conda:
        "../envs/Kraken2.yaml"
    localrule: True
    shell:
        """
        mkdir -p $(dirname {output.versions})
        printf "kraken2\\t%s\\n" "$(kraken2 --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+(\\.[0-9]+)?' | head -1 || echo unknown)" > {output.versions}
        """
