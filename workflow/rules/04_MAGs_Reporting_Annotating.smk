# MAGs REPORTING AND ANNOTATION SNAKEFILE
# Last Updated: 2026-06-06
# Author: Matias Becker Burgos (Matias.BeckerBurgos@unil.ch)

#######################################################
# DESCRIPTION
#######################################################
# MAGs Reporting and annotating rules usable in many different pipelines to classify MAGs created.
#
# MODULE INTERFACE
#   INPUT:       HIGHQC_BINS_DIR (from 03_binning.smk or 03b_Binning_MetaBAT2.smk,
#                depending on the alternative_binning steps flag — set once in Snakefile)
#   OUTPUT:      {OUTPUT_DIR}/02_Results/04_MAGs_Reporting_Annotating/MAG_Representative_info_parsed.tsv
#                {OUTPUT_DIR}/02_Results/04_MAGs_Reporting_Annotating/reference_SGBs/  (used by 05a)
#   REQUIRES:    03_binning.smk or 03b_Binning_MetaBAT2.smk (HIGHQC_BINS_DIR directory)
#   STEPS KEY:   MAGs_Reporting_step: true
#   USES:        CONTIGS (blank samples excluded — same set as binning)
#   NOTE:        Module 05a (InStrain) depends on this module's dRep outputs.

# Steps included:
# 01: Taxonomic Classification (GTDB-Tk)
# 02: MAGs Tree inference (GTDB-Tk)
# 03: iTOL Parsing
# 04: Dereplication (dRep)
# 05: Reporting Summary


##################################################
# RULES
##################################################

########## 01 TAXONOMIC CLASSIFICATION ###########

rule gtdb_classify:
    input:
        db= config["GTDB_Tk"]["db"],
        filtered_mags= HIGHQC_BINS_DIR
    output:
        # classify_wf's own output tree (identify/align/classify + pplacer scratch)
        # is large and only two files inside it are ever needed again — those are
        # copied out to gtdb_classify_summary/ (below) so this bulk tree can be
        # temp()'d without losing them.
        classified_out= temp(directory(f"{OUTPUT_DIR}/01_Analysis/04_MAGs_Reporting_Annotating/gtdb_classify/")),
        GTDB_Tk= f"{OUTPUT_DIR}/01_Analysis/04_MAGs_Reporting_Annotating/gtdb_classify_summary/gtdbtk.bac120.summary.tsv",
        msa = temp(f"{OUTPUT_DIR}/01_Analysis/04_MAGs_Reporting_Annotating/gtdb_classify_summary/gtdbtk.bac120.user_msa.fasta.gz")
    log:
        f"{LOG_DIR}/04_MAGs_Reporting_Annotating/GTDB_Tk/GTDB_classify.log"
    benchmark:
        f"{BENCH_DIR}/04_MAGs_Reporting_Annotating/GTDB_Tk/GTDB_classify.tsv"
    conda:
        "../envs/Reporting.yaml"
    threads: config["GTDB_Tk"]["threads"]
    params:
        tmp = f"{SCRATCH_DIR}/gtdbtk_classify/",
        scratch = f"{SCRATCH_DIR}/gtdbtk_classify_test_pplacer"
    resources:
        mem_mb= config["GTDB_Tk"]["memory_mb"],
        cpus_per_task= config["GTDB_Tk"]["threads"],
        runtime= config["GTDB_Tk"]["runtime_min"]
    shell:
        """
        export GTDBTK_DATA_PATH={input.db}
        echo "GTDB_Tk DB is in: $GTDBTK_DATA_PATH" >> {log}

        mkdir -p {params.tmp}
        
        gtdbtk classify_wf --genome_dir {input.filtered_mags} \
        --extension fa \
        --tmpdir {params.tmp} \
        --pplacer_cpus 4 \
        --scratch_dir {params.scratch} \
        --out_dir {output.classified_out} \
        --cpus {threads} > {log} 2>&1

        # Promote the two files other rules actually need out of the bulk
        # classify_wf tree before it gets cleaned up (temp()).
        mkdir -p $(dirname {output.GTDB_Tk})
        cp {output.classified_out}/classify/gtdbtk.bac120.summary.tsv {output.GTDB_Tk}
        cp {output.classified_out}/align/gtdbtk.bac120.user_msa.fasta.gz {output.msa}
        """

################ 02 MAGs TREE INFERENCE ###################

rule gtdb_infer:
    input:
        msa = f"{OUTPUT_DIR}/01_Analysis/04_MAGs_Reporting_Annotating/gtdb_classify_summary/gtdbtk.bac120.user_msa.fasta.gz",
        db= config["GTDB_Tk"]["db"]
    output:
        # Only gtdb_decorate reads one file out of this (IQ-TREE working dir);
        # nothing needs the rest once that's extracted.
        out_dir = temp(directory(f"{OUTPUT_DIR}/01_Analysis/04_MAGs_Reporting_Annotating/gtdb_infer/bac120_infer_tree"))
    log:
        f"{LOG_DIR}/04_MAGs_Reporting_Annotating/GTDB_Tk/GTDB_infer.log"
    benchmark:
        f"{BENCH_DIR}/04_MAGs_Reporting_Annotating/GTDB_Tk/GTDB_infer.tsv"
    conda:
        "../envs/Reporting.yaml"
    threads: config["GTDB_Tk"]["threads"]
    resources:
        mem_mb= config["GTDB_Tk"]["memory_mb"],
        cpus_per_task= config["GTDB_Tk"]["threads"],
        runtime= config["GTDB_Tk"]["runtime_min"]
    shell:
        """
        export GTDBTK_DATA_PATH={input.db}
        echo "GTDB_Tk DB is in: $GTDBTK_DATA_PATH" >> {log}

        gtdbtk infer --msa_file {input.msa} \
        --cpus {threads} \
        --out_dir {output.out_dir} > {log} 2>&1
        """

rule gtdb_decorate:
    input:
        tree_dir= f"{OUTPUT_DIR}/01_Analysis/04_MAGs_Reporting_Annotating/gtdb_infer/bac120_infer_tree",
        gtdb_classify= f"{OUTPUT_DIR}/01_Analysis/04_MAGs_Reporting_Annotating/gtdb_classify_summary/gtdbtk.bac120.summary.tsv",
        db= config["GTDB_Tk"]["db"]
    output:
        tree= temp(f"{OUTPUT_DIR}/01_Analysis/04_MAGs_Reporting_Annotating/gtdb_decorate/gtdbtk.unrooted.decorated.tree")
    log:
        f"{LOG_DIR}/04_MAGs_Reporting_Annotating/GTDB_Tk/GTDB_decorate.log"
    benchmark:
        f"{BENCH_DIR}/04_MAGs_Reporting_Annotating/GTDB_Tk/GTDB_decorate.tsv"
    conda:
        "../envs/Reporting.yaml"
    threads: 1
    localrule: True
    shell:
        """
        export GTDBTK_DATA_PATH={input.db}
        echo "GTDB_Tk DB is in: $GTDBTK_DATA_PATH" >> {log}

        gtdbtk decorate --input_tree {input.tree_dir}/gtdbtk.unrooted.tree \
        --output_tree {output.tree} \
        --gtdbtk_classification_file {input.gtdb_classify} > {log} 2>&1
        """

################ 03 iTOL PARSING #########################

rule convert_to_iTOL:
    input:
        db = config["GTDB_Tk"]["db"],
        tree = f"{OUTPUT_DIR}/01_Analysis/04_MAGs_Reporting_Annotating/gtdb_decorate/gtdbtk.unrooted.decorated.tree"
    output:
        iTOL = f"{OUTPUT_DIR}/02_Results/04_MAGs_Reporting_Annotating/iTOL_visualisation/gtdbtk.bac120.iTOL.MAGs.tree"
    log:
        f"{LOG_DIR}/04_MAGs_Reporting_Annotating/iTOL_visualisation/iTOL_tree_formating.log"
    benchmark:
        f"{BENCH_DIR}/04_MAGs_Reporting_Annotating/iTOL_visualisation/iTOL_tree_formating.tsv"
    conda:
        "../envs/Reporting.yaml"
    threads: 1
    localrule: True
    shell:
        """
        export GTDBTK_DATA_PATH={input.db}
        echo "GTDB_Tk DB is in: $GTDBTK_DATA_PATH" >> {log}

        gtdbtk convert_to_itol --input_tree {input.tree} --output_tree {output.iTOL} > {log} 2>&1
        """


rule parse_gtdbtk_summary_to_iTOL:
    input:
        gtdb_summary = f"{OUTPUT_DIR}/01_Analysis/04_MAGs_Reporting_Annotating/gtdb_classify_summary/gtdbtk.bac120.summary.tsv"
    output:
        dir = directory(f"{OUTPUT_DIR}/02_Results/04_MAGs_Reporting_Annotating/iTOL_visualisation/Annotation_files")
    params:
        metadata= config["DATA_VARIABLES"]["metadata_path"],
        config= "../config/config.yaml"
    log:
        f"{LOG_DIR}/04_MAGs_Reporting_Annotating/iTOL_visualisation/iTOL_Annotation_files.log"
    benchmark:
        f"{BENCH_DIR}/04_MAGs_Reporting_Annotating/iTOL_visualisation/iTOL_Annotation_files.tsv"
    conda:
        "../envs/Parsing.yaml"
    threads: 1
    localrule: True
    shell:
        """
        python3 scripts/parse_gtdbtk_to_itol.py --metadata {params.metadata} \
        --summary {input.gtdb_summary} \
        --config {params.config} \
        --outdir {output.dir} > {log} 2>&1
        """


rule filter_gtdbtk_for_multiqc:
    input:
        gtdb_summary = f"{OUTPUT_DIR}/01_Analysis/04_MAGs_Reporting_Annotating/gtdb_classify_summary/gtdbtk.bac120.summary.tsv"
    output:
        filtered = f"{OUTPUT_DIR}/02_Results/04_MAGs_Reporting_Annotating/gtdbtk_summary_multiqc.tsv"
    conda:
        "../envs/Parsing.yaml"
    log:
        f"{LOG_DIR}/04_MAGs_Reporting_Annotating/filter_gtdbtk_for_multiqc.log"
    benchmark:
        f"{BENCH_DIR}/04_MAGs_Reporting_Annotating/filter_gtdbtk_for_multiqc.tsv"
    threads: 1
    localrule: True
    shell:
        """
        python3 scripts/filter_gtdbtk_summary_for_multiqc.py \
        --input {input.gtdb_summary} \
        --output {output.filtered} > {log} 2>&1
        """

#################### 04 DEREPLICATION #############################

rule parse_dRep_genomeInfo:
    input:
        CheckM_summary = f"{OUTPUT_DIR}/02_Results/03_Binning/CheckM_QC/filtered_MAGs_stats.tsv"
    output:
        GenomeInfo = temp(f"{OUTPUT_DIR}/02_Results/04_MAGs_Reporting_Annotating/GenomeInfo_dRep.csv")
    conda:
        "../envs/Parsing.yaml"
    log:
        f"{LOG_DIR}/04_MAGs_Reporting_Annotating/parse_dRep_genomeInfo/parse_dRep_genomeInfo.log"
    benchmark:
        f"{BENCH_DIR}/04_MAGs_Reporting_Annotating/parse_dRep_genomeInfo/parse_dRep_genomeInfo.tsv"
    threads: 1
    localrule: True
    shell:
        """
        python3 scripts/parse_dRep_genomeInfo.py \
        --CheckM_summary {input.CheckM_summary} \
        --output {output.GenomeInfo}
        """


rule run_dRep:
    input:
        mag = HIGHQC_BINS_DIR,
        GenomeInfo = f"{OUTPUT_DIR}/02_Results/04_MAGs_Reporting_Annotating/GenomeInfo_dRep.csv"
    output:
        drep_out= directory(f"{OUTPUT_DIR}/02_Results/04_MAGs_Reporting_Annotating/reference_SGBs/")
    conda:
        "../envs/Reporting.yaml"
    params:
        completeness= config["dRep"]["completeness_threshold"],
        contamination= config["dRep"]["contamination_threshold"],
        ANI= config["dRep"]["secondary_clustering_ANI"]
    log:
        f"{LOG_DIR}/04_MAGs_Reporting_Annotating/dRep/dRep_report.log"
    benchmark:
        f"{BENCH_DIR}/04_MAGs_Reporting_Annotating/dRep/dRep_report.tsv"
    threads: config["dRep"]["threads"]
    resources:
        mem_mb = config["dRep"]["memory_mb"],
        cpus_per_task = config["dRep"]["threads"],
        runtime = config["dRep"]["runtime_min"]
    shell:
        """
        dRep dereplicate {output.drep_out} \
            -g {input.mag}/* \
            --genomeInfo {input.GenomeInfo} \
            --completeness {params.completeness} \
            --contamination {params.contamination} \
            -sa {params.ANI} \
            -p {threads} \
            -d > {log} 2>&1
        """

########################## 05 REPORTING SUMMARY ##############################

rule MAG_Representative_info:
    input:
        CheckM= f"{OUTPUT_DIR}/02_Results/03_Binning/CheckM_QC/filtered_MAGs_stats.tsv",
        GTDB_Tk= f"{OUTPUT_DIR}/01_Analysis/04_MAGs_Reporting_Annotating/gtdb_classify_summary/gtdbtk.bac120.summary.tsv",
        dRep= f"{OUTPUT_DIR}/02_Results/04_MAGs_Reporting_Annotating/reference_SGBs/"
    output:
        Summary= f"{OUTPUT_DIR}/02_Results/04_MAGs_Reporting_Annotating/MAG_Representative_info_parsed.tsv"
    conda:
        "../envs/Parsing.yaml"
    log:
        f"{LOG_DIR}/04_MAGs_Reporting_Annotating/MAG_Representative_info.log"
    benchmark:
        f"{BENCH_DIR}/04_MAGs_Reporting_Annotating/MAG_Representative_info.tsv"
    threads: 1
    localrule: True
    shell:
        """
        python3 scripts/parse_MAG_representative_info.py \
        --CheckM_summary {input.CheckM} \
        --GTDB_Tk_summary {input.GTDB_Tk} \
        --drep_dir {input.dRep} \
        --output {output.Summary} > {log} 2>&1
        """