# BINNING SNAKEFILE
# Last Updated: 2026-01-14
# Author: Matias Becker Burgos (Matias.BeckerBurgos@unil.ch)

#######################################################
# DESCRIPTION
#######################################################
# Binning rules usable in many different pipelines
# The workflow start with the assumption that all the data files are in the same directory.

# Steps included:
# 01: Fairy sketching
# 02: Fairy coverage
# 03: Mixte Binning
#   03a: Metabat2
#   03b: Maxbin2
#   03c: SemiBin2
# 04: Bin Selection (Das Tool)
# 05: Bin QC (CheckM2)
# 06: Bin QC filtering


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

#################### 02 FAIRY COVERAGE #########################


#################### 03A METABAT2 ####################

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
        Bins = temp(directory(f"{OUTPUT_DIR}/01_Analysis/03_Binning/MetaBat2/Bins/{{contig}}_metabat2"))
    params:
        min_contig_length = config["MetaBat2"]["min_contig_length"],
        basename = lambda wc: f"{OUTPUT_DIR}/01_Analysis/03_Binning/MetaBat2/Bins/{wc.contig}_metabat2/{wc.contig}_metabat2.bin"
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
                 --minContig {params.min_contig_length} 2> {log} 1>&2
        """


######################### 03B MAXBIN2 ################################

rule Fairy_coverage_maxbin2:
    input:
        sketches = expand(f"{SCRATCH_DIR}/Fairy_Sketch/{{sample}}.paired.bcsp", sample=CONTIGS),
        contig = (
            f"{OUTPUT_DIR}/02_Results/02_Assembly/Contig_Filtering/"
            f"{{contig}}_min{config['Contig_filter']['length_threshold']}_contigs.fasta"
        )
    output:
        wide = temp(f"{OUTPUT_DIR}/01_Analysis/03_Binning/MaxBin2/Fairy_coverage/{{contig}}_maxbin_wide.tsv")
    conda:
        "../envs/Fairy.yaml"
    threads: config["Fairy_coverage"]["threads"]
    log:
        f"{LOG_DIR}/03_Binning/Fairy/{{contig}}_maxbin_coverage.log"
    benchmark:
        f"{BENCH_DIR}/03_Binning/Fairy/{{contig}}_maxbin_coverage.tsv"
    resources:
        mem_mb = config["Fairy_coverage"]["memory_mb"],
        runtime = config["Fairy_coverage"]["runtime_min"],
        cpus_per_task = config["Fairy_coverage"]["threads"]
    shell:
        """
        fairy coverage {input.sketches} {input.contig} \
            --maxbin-format \
            -t {threads} \
            -o {output.wide} > {log} 2>&1
        """


rule maxbin2_abund_list:
    input:
        wide = f"{OUTPUT_DIR}/01_Analysis/03_Binning/MaxBin2/Fairy_coverage/{{contig}}_maxbin_wide.tsv"
    output:
        abund_list = temp(f"{OUTPUT_DIR}/01_Analysis/03_Binning/MaxBin2/Abundance_list/{{contig}}_abund_list.txt"),
        per_sample = temp(expand(
            f"{SCRATCH_DIR}/Coverage_table/{{{{contig}}}}/{{sample}}_coverage.txt", sample=CONTIGS
        ))
    params:
        outdir = lambda wc: f"{SCRATCH_DIR}/Coverage_table/{wc.contig}",
        samples = " ".join(CONTIGS),
    conda:
        "../envs/Parsing.yaml"
    benchmark:
        f"{BENCH_DIR}/03_Binning/MaxBin2/Abundance_list/{{contig}}_abund_list.tsv"
    log:
        f"{LOG_DIR}/03_Binning/MaxBin2/Abundance_list/{{contig}}_abund_list.log"
    shell:
        """
        mkdir -p {params.outdir}
        python scripts/split_fairy_maxbin_table.py \
            --input {input.wide} \
            --outdir {params.outdir} \
            --samples {params.samples} \
            --output-list {output.abund_list} > {log} 2>&1
        """


rule MaxBin2:
    input:
        contigs = (
            f"{OUTPUT_DIR}/02_Results/02_Assembly/Contig_Filtering/"
            f"{{contig}}_min{config['Contig_filter']['length_threshold']}_contigs.fasta"
        ),
        abund_list = f"{OUTPUT_DIR}/01_Analysis/03_Binning/MaxBin2/Abundance_list/{{contig}}_abund_list.txt"
    output:
        bins = temp(directory(f"{OUTPUT_DIR}/01_Analysis/03_Binning/MaxBin2/Bins/{{contig}}_maxbin2"))
    params:
        basename = lambda wc: f"{OUTPUT_DIR}/01_Analysis/03_Binning/MaxBin2/Bins/{wc.contig}_maxbin2/{wc.contig}_maxbin2.bin",
        min_contig_length = config["MaxBin2"]["min_contig_length"]
    threads: config["MaxBin2"]["threads"]
    conda:
        "../envs/Binning.yaml"
    benchmark:
        f"{BENCH_DIR}/03_Binning/MaxBin2/Binning/{{contig}}_maxbin2.log"
    log:
        f"{LOG_DIR}/03_Binning/MaxBin2/Binning/{{contig}}_maxbin2.log"
    resources:
        mem_mb = config["MaxBin2"]["memory_mb"],
        runtime = config["MaxBin2"]["runtime_min"],
        cpus_per_task = config["MaxBin2"]["threads"]
    shell:
        """
        mkdir -p {output.bins}

        run_MaxBin.pl -thread {threads} \
        -min_contig_length {params.min_contig_length} \
        -contig {input.contigs} \
        -abund_list {input.abund_list} \
        -out {params.basename} 2> {log} 1>&2

        # Normalize MaxBin2 bin names
        for f in {output.bins}/*.fasta; do
            bn=$(basename "$f" .fasta)
            id=$(echo "$bn" | awk -F. '{{printf "%d", $NF}}')
            mv "$f" "{output.bins}/{{wc.contig}}_maxbin2.bin.${{id}}.fa"
        done

        """


######################### 03C SEMIBIN2 ################################

rule Fairy_coverage_aemb:
    input:
        sketches = expand(f"{SCRATCH_DIR}/Fairy_Sketch/{{sample}}.paired.bcsp", sample=CONTIGS),
        contig = (
            f"{OUTPUT_DIR}/02_Results/02_Assembly/Contig_Filtering/"
            f"{{contig}}_min{config['Contig_filter']['length_threshold']}_contigs.fasta"
        )
    output:
        # One aemb file per sample: fairy coverage --aemb-format only produces a true
        # single-column-per-sample file when given one sketch at a time (confirmed —
        # given all sketches at once it instead produces the same wide table as
        # --maxbin-format, just without a header). Each call is near-instant (no
        # alignment), so N per-contig calls stay cheap despite the O(N) loop.
        aemb_dir = temp(directory(f"{OUTPUT_DIR}/01_Analysis/03_Binning/SemiBin2/Fairy_aemb/{{contig}}"))
    params:
        sketch_dir = f"{SCRATCH_DIR}/Fairy_Sketch",
        samples = CONTIGS
    conda:
        "../envs/Fairy.yaml"
    threads: config["Fairy_coverage"]["threads"]
    log:
        f"{LOG_DIR}/03_Binning/Fairy/{{contig}}_aemb_coverage.log"
    benchmark:
        f"{BENCH_DIR}/03_Binning/Fairy/{{contig}}_aemb_coverage.tsv"
    resources:
        mem_mb = config["Fairy_coverage"]["memory_mb"],
        runtime = config["Fairy_coverage"]["runtime_min"],
        cpus_per_task = config["Fairy_coverage"]["threads"]
    shell:
        """
        mkdir -p {output.aemb_dir}
        for sample in {params.samples}; do
            fairy coverage {params.sketch_dir}/${{sample}}.paired.bcsp {input.contig} \
                --aemb-format -t {threads} \
                -o {output.aemb_dir}/${{sample}}_aemb.txt
        done > {log} 2>&1
        """


rule SemiBin2_binning:
    input:
        contig = (
            f"{OUTPUT_DIR}/02_Results/02_Assembly/Contig_Filtering/"
            f"{{contig}}_min{config['Contig_filter']['length_threshold']}_contigs.fasta"
        ),
        aemb_dir = f"{OUTPUT_DIR}/01_Analysis/03_Binning/SemiBin2/Fairy_aemb/{{contig}}"
    output:
        bins = temp(directory(f"{OUTPUT_DIR}/01_Analysis/03_Binning/SemiBin2/Bins/{{contig}}_semibin2"))
    conda:
        "../envs/SemiBin2.yaml"
    threads: config["SemiBin2"]["threads"]
    log:
        f"{LOG_DIR}/03_Binning/SemiBin2/{{contig}}_semibin2.log"
    benchmark:
        f"{BENCH_DIR}/03_Binning/SemiBin2/{{contig}}_semibin2.tsv"
    resources:
        mem_mb = config["SemiBin2"]["memory_mb"],
        runtime = config["SemiBin2"]["runtime_min"],
        cpus_per_task = config["SemiBin2"]["threads"]
    shell:
        """
        SemiBin2 single_easy_bin \
            -i {input.contig} \
            -a {input.aemb_dir}/*.txt \
            -o {output.bins} \
            --compression none \
            -t {threads} > {log} 2>&1
        """


rule tool_version_semibin2:
    output:
        versions = temp(f"{OUTPUT_DIR}/01_Analysis/ToolVersions/semibin2_versions.tsv")
    conda:
        "../envs/SemiBin2.yaml"
    localrule: True
    shell:
        """
        mkdir -p $(dirname {output.versions})
        printf "semibin2\\t%s\\n" "$(SemiBin2 --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+(\\.[0-9]+)?' | head -1 || echo unknown)" > {output.versions}
        """

######################### 04 BIN SELECTION ###########################

rule metabat2_Scaffolds2Bin:
    input:
        bins = f"{OUTPUT_DIR}/01_Analysis/03_Binning/MetaBat2/Bins/{{contig}}_metabat2/"
    output:
        scaffolds2bin = temp(f"{OUTPUT_DIR}/01_Analysis/03_Binning/MetaBat2/scaffolds2bin/{{contig}}_Scaffolds2Bin.tsv")
    conda:
        "../envs/Binning.yaml"
    benchmark:
        f"{BENCH_DIR}/03_Binning/MetaBat2/scaffolds2bin/{{contig}}_Scaffolds2Bin.tsv"
    log:
        f"{LOG_DIR}/03_Binning/MetaBat2/scaffolds2bin/{{contig}}_Scaffolds2Bin.log"
    shell:
        """
            Fasta_to_Contig2Bin.sh \
            -i {input.bins} \
            -e fa \
            | awk 'BEGIN{{OFS="\t"}} {{print $1,$NF}}' > {output.scaffolds2bin}
        """


rule maxbin2_Scaffolds2Bin:
    input:
        bins = f"{OUTPUT_DIR}/01_Analysis/03_Binning/MaxBin2/Bins/{{contig}}_maxbin2/"
    output:
        scaffolds2bin = temp(f"{OUTPUT_DIR}/01_Analysis/03_Binning/MaxBin2/scaffolds2bin/{{contig}}_Scaffolds2Bin.tsv")
    conda:
        "../envs/Binning.yaml"
    benchmark:
        f"{BENCH_DIR}/03_Binning/MaxBin2/scaffolds2bin/{{contig}}_Scaffolds2Bin.tsv"
    log:
        f"{LOG_DIR}/03_Binning/MaxBin2/scaffolds2bin/{{contig}}_Scaffolds2Bin.log"
    shell:
        """
        Fasta_to_Contig2Bin.sh \
        -i {input.bins} \
        -e fasta > {output.scaffolds2bin}
        """


rule semibin2_Scaffolds2Bin:
    input:
        bins = f"{OUTPUT_DIR}/01_Analysis/03_Binning/SemiBin2/Bins/{{contig}}_semibin2"
    output:
        scaffolds2bin = temp(f"{OUTPUT_DIR}/01_Analysis/03_Binning/SemiBin2/scaffolds2bin/{{contig}}_Scaffolds2Bin.tsv")
    conda:
        "../envs/Binning.yaml"
    benchmark:
        f"{BENCH_DIR}/03_Binning/SemiBin2/scaffolds2bin/{{contig}}_Scaffolds2Bin.tsv"
    log:
        f"{LOG_DIR}/03_Binning/SemiBin2/scaffolds2bin/{{contig}}_Scaffolds2Bin.log"
    shell:
        """
        # NOTE: output_bins/ is SemiBin2's documented default results subdirectory —
        # not independently confirmed in this session (our synthetic test fixture has
        # too few contigs for SemiBin2's clustering to complete past input validation).
        # Verify this path on the first real run against actual assembled contigs.
        Fasta_to_Contig2Bin.sh \
        -i {input.bins}/output_bins \
        -e fa \
        | awk 'BEGIN{{OFS="\t"}} {{print $1,$NF}}' > {output.scaffolds2bin} 2> {log}
        """


rule DAS_Tool:
    input:
        metabat2 = f"{OUTPUT_DIR}/01_Analysis/03_Binning/MetaBat2/scaffolds2bin/{{contig}}_Scaffolds2Bin.tsv",
        maxbin2 = f"{OUTPUT_DIR}/01_Analysis/03_Binning/MaxBin2/scaffolds2bin/{{contig}}_Scaffolds2Bin.tsv",
        semibin2 = f"{OUTPUT_DIR}/01_Analysis/03_Binning/SemiBin2/scaffolds2bin/{{contig}}_Scaffolds2Bin.tsv",
        contigs = (
            f"{OUTPUT_DIR}/02_Results/02_Assembly/Contig_Filtering/"
            f"{{contig}}_min{config['Contig_filter']['length_threshold']}_contigs.fasta"
        )
    output:
        out = f"{OUTPUT_DIR}/01_Analysis/03_Binning/DAS_Tool/{{contig}}_DASTool/{{contig}}_DASTool_summary.tsv",
        dir = directory(f"{OUTPUT_DIR}/01_Analysis/03_Binning/DAS_Tool/{{contig}}_DASTool/{{contig}}_DASTool_bins")
    params:
        basename = lambda wc: f"{OUTPUT_DIR}/01_Analysis/03_Binning/DAS_Tool/{wc.contig}_DASTool/{wc.contig}",
        search_engine = config["DASTool"]["search_engine"]
    conda:
        "../envs/Binning.yaml"
    threads: config["DASTool"]["threads"]
    benchmark:
        f"{BENCH_DIR}/03_Binning/Das_Tool/{{contig}}_binning.tsv"
    log:
        f"{LOG_DIR}/03_Binning/Das_Tool/{{contig}}_binning.log"
    resources:
        mem_mb = config["DASTool"]["memory_mb"],
        cpus_per_task = config["DASTool"]["threads"],
        runtime = config["DASTool"]["runtime_min"]
    shell:
        """
        DAS_Tool \
        --bins {input.metabat2},{input.maxbin2},{input.semibin2} \
        --contigs {input.contigs} \
        --outputbasename {params.basename} \
        --labels metabat2,maxbin2,semibin2 \
        --write_bins \
        --write_bin_evals \
        --threads {threads} \
        --search_engine {params.search_engine}

        # Rename all bins to consistent pattern DAS_Tool as inconsistent bin naming
        i=1
        for f in {output.dir}/*.fa; do
            mv -n "$f" "{output.dir}/{wildcards.contig}_DASTool_bin_${{i}}.fa"
            i=$((i+1))
        done
        """


rule aggregate_dastool_summary:
    input:
        summaries = expand(
            f"{OUTPUT_DIR}/01_Analysis/03_Binning/DAS_Tool/{{contig}}_DASTool/{{contig}}_DASTool_summary.tsv",
            contig=CONTIGS,
        )
    output:
        summary = f"{OUTPUT_DIR}/02_Results/03_Binning/DASTool_binning_summary.tsv"
    params:
        samples = CONTIGS
    conda:
        "../envs/Parsing.yaml"
    log:
        f"{LOG_DIR}/03_Binning/Das_Tool/aggregate_dastool_summary.log"
    benchmark:
        f"{BENCH_DIR}/03_Binning/Das_Tool/aggregate_dastool_summary.tsv"
    threads: 1
    localrule: True
    shell:
        """
        python scripts/aggregate_dastool_summary.py \
        --input {input.summaries} \
        --samples {params.samples} \
        --output {output.summary} > {log} 2>&1
        """


######################### 05 BIN QC ##################################

rule CheckM_QC:
    input:
        summary = f"{OUTPUT_DIR}/01_Analysis/03_Binning/DAS_Tool/{{contig}}_DASTool/{{contig}}_DASTool_summary.tsv",
        dir = f"{OUTPUT_DIR}/01_Analysis/03_Binning/DAS_Tool/{{contig}}_DASTool/{{contig}}_DASTool_bins"
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

######################### 06 BIN QC FILTERING ##################################

rule Bin_filtering_QC_summary:
    input:
        stats = expand(f"{OUTPUT_DIR}/01_Analysis/03_Binning/CheckM_QC/{{contig}}_checkm_QC/quality_report.tsv", contig=CONTIGS),
        bins = expand(f"{OUTPUT_DIR}/01_Analysis/03_Binning/DAS_Tool/{{contig}}_DASTool/{{contig}}_DASTool_bins", contig=CONTIGS)
    output:
        full_stats= f"{OUTPUT_DIR}/02_Results/03_Binning/CheckM_QC/all_MAGs_stats.tsv",
        filtered_stats= f"{OUTPUT_DIR}/02_Results/03_Binning/CheckM_QC/filtered_MAGs_stats.tsv",
        filtered_mags= directory(f"{OUTPUT_DIR}/02_Results/03_Binning/HighQC_Bins/")
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
        -e {params.contamination} >> {log}
        """
