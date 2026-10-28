# BINNING METABAT2 SNAKEFILE (ALTERNATIVE STANDALONE)
# Last Updated: 2026-06-06
# Author: Matias Becker Burgos (Matias.BeckerBurgos@unil.ch)

#######################################################
# DESCRIPTION
#######################################################
# Standalone MetaBAT2 binning — alternative to the full consensus binning in 03_binning.smk.
# Run this OR 03_binning.smk, not both simultaneously (outputs go to different directories).
#
# MODULE INTERFACE
#   INPUT:       Filtered contigs from 02_assembly.smk
#   OUTPUT:      {OUTPUT_DIR}/02_Results/03_Binning/HighQC_Bins_MetaBAT2/
#   REQUIRES:    02_assembly.smk (Contig_Filtering)
#   STEPS KEY:   alternative_binning: true
#   USES:        CONTIGS (blank samples excluded)
#   NOTE:        Writes to HighQC_Bins_MetaBAT2/ (NOT HighQC_Bins/) to avoid
#                conflicting with 03_binning.smk output. Module 04 reads from
#                HIGHQC_BINS_DIR (set once in Snakefile from this same
#                alternative_binning flag), so no separate config update is needed.

# Steps included:
# 01: Fairy sketching
# 02: Fairy coverage
# 03: Metabat2
# 04: Bin QC (CheckM2)
# 05: Bin QC filtering


##################################################
# RULES
##################################################

#################### 01 FAIRY SKETCHING ####################

rule Fairy_sketch:
    input:
        R1 = expand(f"{CLEAN_DATA_DIR}/{{sample}}_R1_cleaned.fastq.gz", sample=CONTIGS),
        R2 = expand(f"{CLEAN_DATA_DIR}/{{sample}}_R2_cleaned.fastq.gz", sample=CONTIGS),
    output:
        sketches = temp(expand(f"{SCRATCH_DIR}/Fairy_Sketch/{{sample}}.paired.bcsp", sample=CONTIGS))
    params:
        sketch_dir = f"{SCRATCH_DIR}/Fairy_Sketch",
        samples = " ".join(CONTIGS),
    conda:
        "../envs/Fairy.yaml"
    threads: config["Fairy_sketch"]["threads"]
    log:
        f"{LOG_DIR}/03_Binning/Fairy/sketch.log"
    benchmark:
        f"{BENCH_DIR}/03_Binning/Fairy/sketch.tsv"
    resources:
        mem_mb = config["Fairy_sketch"]["memory_mb"],
        runtime = config["Fairy_sketch"]["runtime_min"],
        cpus_per_task = config["Fairy_sketch"]["threads"]
    shell:
        """
        mkdir -p {params.sketch_dir}
        fairy sketch -1 {input.R1} -2 {input.R2} \
            -S {params.samples} \
            -d {params.sketch_dir} \
            -t {threads} > {log} 2>&1
        """


rule tool_version_fairy:
    output:
        versions = temp(f"{OUTPUT_DIR}/01_Analysis/ToolVersions/fairy_versions.tsv")
    conda:
        "../envs/Fairy.yaml"
    localrule: True
    shell:
        """
        mkdir -p $(dirname {output.versions})
        printf "fairy\\t%s\\n" "$(fairy --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+(\\.[0-9]+)?' | head -1 || echo unknown)" > {output.versions}
        """

#################### 02 FAIRY COVERAGE ####################

rule Fairy_coverage_metabat2:
    input:
        sketches = expand(f"{SCRATCH_DIR}/Fairy_Sketch/{{sample}}.paired.bcsp", sample=CONTIGS),
        contig = (
            f"{OUTPUT_DIR}/02_Results/02_Assembly/Contig_Filtering/"
            f"{{contig}}_min{config['Contig_filter']['length_threshold']}_contigs.fasta"
        )
    output:
        depth = temp(f"{OUTPUT_DIR}/01_Analysis/03_Binning/MetaBat2/Depth_Contig/{{contig}}_depth.txt")
    conda:
        "../envs/Fairy.yaml"
    threads: config["Fairy_coverage"]["threads"]
    log:
        f"{LOG_DIR}/03_Binning/Fairy/{{contig}}_coverage.log"
    benchmark:
        f"{BENCH_DIR}/03_Binning/Fairy/{{contig}}_coverage.tsv"
    resources:
        mem_mb = config["Fairy_coverage"]["memory_mb"],
        runtime = config["Fairy_coverage"]["runtime_min"],
        cpus_per_task = config["Fairy_coverage"]["threads"]
    shell:
        """
        fairy coverage {input.sketches} {input.contig} \
            -t {threads} \
            -o {output.depth} > {log} 2>&1
        """


rule Metabat2:
    input:
        contig = (
            f"{OUTPUT_DIR}/02_Results/02_Assembly/Contig_Filtering/"
            f"{{contig}}_min{config['Contig_filter']['length_threshold']}_contigs.fasta"
        ),
        depth = f"{OUTPUT_DIR}/01_Analysis/03_Binning/MetaBat2/Depth_Contig/{{contig}}_depth.txt"
    output:
        Bins = directory(f"{OUTPUT_DIR}/01_Analysis/03_Binning/MetaBat2/Bins/{{contig}}_metabat2")
    params:
        min_contig_length = config["MetaBat2"]["min_contig_length"],
        basename = lambda wc: f"{OUTPUT_DIR}/01_Analysis/03_Binning/MetaBat2/Bins/{wc.contig}_metabat2/{wc.contig}_bin"
    conda:
        "../envs/Binning.yaml"
    threads: config["MetaBat2"]["threads"]
    log:
        f"{LOG_DIR}/03_Binning/MetaBat2/Binning/{{contig}}_metabat2.log"
    benchmark:
        f"{BENCH_DIR}/03_Binning/MetaBat2/Binning/{{contig}}_metabat2.tsv"
    resources:
        mem_mb = config["MetaBat2"]["memory_mb"],
        runtime = config["MetaBat2"]["runtime_min"],
        cpus_per_task = config["MetaBat2"]["threads"]
    shell:
        """
        mkdir -p {output.Bins}
        metabat2 --numThreads {threads} \
                 --inFile {input.contig} \
                 --outFile {params.basename} \
                 --abdFile {input.depth} \
                 --minContig {params.min_contig_length} > {log} 2>&1
        """

######################### 04 BIN QC ##################################

rule CheckM_QC:
    input:
        dir = f"{OUTPUT_DIR}/01_Analysis/03_Binning/MetaBat2/Bins/{{contig}}_metabat2/"
    output:
        file = temp(f"{OUTPUT_DIR}/01_Analysis/03_Binning/CheckM_QC/{{contig}}_checkm_QC/quality_report.tsv")
    params:
        outdir = lambda wc: f"{OUTPUT_DIR}/01_Analysis/03_Binning/CheckM_QC/{wc.contig}_checkm_QC",
        db = config["CheckM"]["db"]
    conda:
        "../envs/Binning.yaml"
    log:
        f"{LOG_DIR}/03_Binning/CheckM/{{contig}}_checkm_QC.log"
    benchmark:
        f"{BENCH_DIR}/03_Binning/CheckM/{{contig}}_checkm_QC.tsv"
    threads: config["CheckM"]["threads"]
    resources:
        mem_mb = config["CheckM"]["memory_mb"],
        runtime = config["CheckM"]["runtime_min"],
        cpus_per_task = config["CheckM"]["threads"]
    shell:
        """
        checkm2 predict --threads {threads} --input {input.dir} --output-directory {params.outdir} \
            -x fa --database_path {params.db}
        """


rule tool_version_binning:
    output:
        versions = temp(f"{OUTPUT_DIR}/01_Analysis/ToolVersions/binning_versions.tsv")
    conda:
        "../envs/Binning.yaml"
    localrule: True
    shell:
        """
        mkdir -p $(dirname {output.versions})
        {{
            printf "metabat2\\t%s\\n" "$(metabat2 --help 2>&1 | grep -oE '[0-9]+\\.[0-9]+(\\.[0-9]+)?' | head -1 || echo unknown)"
            printf "maxbin2\\t%s\\n" "$(run_MaxBin.pl -v 2>&1 | grep -oE '[0-9]+\\.[0-9]+(\\.[0-9]+)?' | head -1 || echo unknown)"
            printf "das_tool\\t%s\\n" "$(DAS_Tool --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+(\\.[0-9]+)?' | head -1 || echo unknown)"
            printf "checkm2\\t%s\\n" "$(checkm2 --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+(\\.[0-9]+)?' | head -1 || echo unknown)"
        }} > {output.versions}
        """

######################### 05 BIN QC FILTERING ##################################

rule Bin_filtering_QC_summary:
    input:
        stats = expand(f"{OUTPUT_DIR}/01_Analysis/03_Binning/CheckM_QC/{{contig}}_checkm_QC/quality_report.tsv", contig=CONTIGS),
        bins = expand(f"{OUTPUT_DIR}/01_Analysis/03_Binning/MetaBat2/Bins/{{contig}}_metabat2", contig=CONTIGS)
    output:
        full_stats= f"{OUTPUT_DIR}/02_Results/03_Binning/CheckM_QC/all_MAGs_stats.tsv",
        filtered_stats= f"{OUTPUT_DIR}/02_Results/03_Binning/CheckM_QC/filtered_MAGs_stats.tsv",
        filtered_mags= directory(f"{OUTPUT_DIR}/02_Results/03_Binning/HighQC_Bins_MetaBAT2/")
    conda:
        "../envs/Parsing.yaml"
    log:
        f"{LOG_DIR}/03_Binning/CheckM/aggregate_checkm.log"
    threads: 1
    localrule: True
    params:
        completeness= config["Bin_filtering"]["completeness"],
        contamination= config["Bin_filtering"]["contamination"]
    shell:
        """
        python scripts/aggregate_checkm.py -i {input.stats} \
        -b {input.bins} \
        -o {output.full_stats} \
        -f {output.filtered_stats} \
        -m {output.filtered_mags} \
        -c {params.completeness} \
        -e {params.contamination}
        """
