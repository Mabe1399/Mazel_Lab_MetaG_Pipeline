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
#   - A Kraken2 database is available at config[kraken2_db][DB_dir] (KRAKEN2_DB)
#
# MODULE INTERFACE
#   INPUT:       Cleaned FASTQ reads in CLEAN_DATA_DIR (from module 01, or pre-cleaned input)
#   OUTPUT:      {OUTPUT_DIR}/02_Results/02b_Taxonomic_Profiling/Kraken2/Combined_kraken2_tax_profiling_report.tsv
#   REQUIRES:    — (works with QC: true or false)
#   STEPS KEY:   taxonomic_profiling: true
#   USES:        SAMPLES (all samples including blanks), KRAKEN2_DB
#   NOTE:        Distinct from module 01's run_kraken2, which classifies the *trimmed*
#                (pre-host-filtering) reads for the final report; this one profiles the
#                cleaned reads. Its scratch .out path must stay different from run_kraken2's,
#                or both rules claim the same file (AmbiguousRuleException with QC: true).

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
        db= KRAKEN2_DB
    output:
        tab=temp(f"{SCRATCH_DIR}/Kraken2_Taxonomic_Profiling/{{sample}}_kraken2.out"),
        rep=temp(f"{OUTPUT_DIR}/01_Analysis/02b_Taxonomic_Profiling/Kraken2/{{sample}}_kraken2.report")
    conda:
        "../envs/Core.yaml"
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
        "../envs/Core.yaml"
    threads: 1
    log:
        f"{LOG_DIR}/02b_Taxonomic_Profiling/Kraken2/Combined_kraken2_parsing.log"
    benchmark:
        f"{BENCH_DIR}/02b_Taxonomic_Profiling/Kraken2/Combined_kraken2_parsing.tsv"
    shell:
        "python3 scripts/Combine_kraken_reports.py -o {output} -i {input} > {log} 2>&1"
