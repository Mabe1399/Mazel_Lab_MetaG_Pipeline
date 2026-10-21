# QUALITY CONTROL AND PREPROCESSING SNAKEFILE
# Last Updated: 2026-06-12
# Author: Matias Becker Burgos (Matias.BeckerBurgos@unil.ch)

#######################################################
# DESCRIPTION
#######################################################
# MODULE INTERFACE
#   INPUT:    Raw FASTQ files in DATA_DIR
#               Multi-lane:       {sample}_L*_R{1,2}*.fastq.gz  (steps: concatenate_lanes: true)
#               Pre-concatenated: {sample}_R{1,2}*.fastq.gz     (steps: concatenate_lanes: false)
#   OUTPUT:   {OUTPUT_DIR}/00_Data/02_Clean_Data/{sample}_R{1,2}_cleaned.fastq.gz
#   REQUIRES: —
#   STEPS KEY: QC
#   USES:     SAMPLES

# Steps included:
# 01: Concatenation (multi-lane) OR symlinking (pre-concatenated)
# 02: Before QC (FastQC + MultiQC)
# 03: Trimming (Trimmomatic)
# 04: After QC (FastQC + MultiQC)
# 05: Taxonomic profiling (Kraken2)
# 06: Host filtering (Bowtie2)

##################################################
# RULES
##################################################

# Derive host index file list from config so it's not hardcoded below
_HOST_INDEX_PREFIX = config["host_filtering"]["index_path"]
_HOST_INDEX_DIR    = os.path.dirname(_HOST_INDEX_PREFIX)
_HOST_INDEX_FILES  = expand(
    f"{_HOST_INDEX_PREFIX}.{{ext}}",
    ext=["1.bt2l", "2.bt2l", "3.bt2l", "4.bt2l", "rev.1.bt2l", "rev.2.bt2l"],
)

_CONCATENATE_LANES = config["steps"].get("concatenate_lanes", True)


def _find_raw_read(wildcards, read):
    """Return the single pre-concatenated FASTQ for a sample/read direction.

    Globs for {sample}_R{read}*.fastq.gz in DATA_DIR.  Raises a clear error
    when no file is found or when multiple files match (ambiguous).
    """
    pattern = f"{DATA_DIR}/{wildcards.sample}_R{read}*.fastq.gz"
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(
            f"concatenate_lanes: false — no R{read} file found for sample "
            f"'{wildcards.sample}' matching {pattern}"
        )
    if len(files) > 1:
        raise ValueError(
            f"concatenate_lanes: false — multiple R{read} files found for "
            f"'{wildcards.sample}': {files}. "
            f"Set concatenate_lanes: true to concatenate them."
        )
    return files[0]


############ 01 Lane handling ######################
# Both branches produce identical output paths so all downstream rules are
# unchanged regardless of the concatenate_lanes setting.

if _CONCATENATE_LANES:
    rule concat_lanes:
        input:
            R1s=lambda wc: sorted(glob.glob(f"{DATA_DIR}/{wc.sample}_L*_R1*.fastq.gz")),
            R2s=lambda wc: sorted(glob.glob(f"{DATA_DIR}/{wc.sample}_L*_R2*.fastq.gz")),
        output:
            # Not temp(): regenerating this after a later invocation finds it
            # missing bumps its mtime, which — since Raw_Trimming/host_filtering
            # read it (directly or transitively) — can cascade into marking
            # already-complete host_filtering/assembly as out of date (assembly
            # is ~23h/sample). See docs/log_Bug_fix.txt.
            R1=f"{SCRATCH_DIR}/Concatenated/{{sample}}_R1_concat.fastq.gz",
            R2=f"{SCRATCH_DIR}/Concatenated/{{sample}}_R2_concat.fastq.gz",
        threads: config["Concatenate"]["threads"]
        log:
            f"{LOG_DIR}/01_QC_Preprocessing/lane_concatenation/{{sample}}_concat.log"
        benchmark:
            f"{BENCH_DIR}/01_QC_Preprocessing/lane_concatenation/{{sample}}_concat.tsv"
        resources:
            runtime = config["Concatenate"]["runtime_min"]
        shell:
            """
            mkdir -p $(dirname {output.R1})
            cat {input.R1s} > {output.R1} 2>  {log}
            cat {input.R2s} > {output.R2} 2>> {log}
            """

else:
    rule link_raw_reads:
        """Symlink pre-concatenated files into the Concatenated/ staging directory."""
        input:
            R1=lambda wc: _find_raw_read(wc, 1),
            R2=lambda wc: _find_raw_read(wc, 2),
        output:
            # Not temp() — see concat_lanes above for why.
            R1=f"{SCRATCH_DIR}/Concatenated/{{sample}}_R1_concat.fastq.gz",
            R2=f"{SCRATCH_DIR}/Concatenated/{{sample}}_R2_concat.fastq.gz",
        localrule: True
        log:
            f"{LOG_DIR}/01_QC_Preprocessing/lane_concatenation/{{sample}}_link.log"
        benchmark:
            f"{BENCH_DIR}/01_QC_Preprocessing/lane_concatenation/{{sample}}_link.tsv"
        shell:
            """
            mkdir -p $(dirname {output.R1})
            ln -sf "$(realpath {input.R1})" {output.R1} > {log} 2>&1
            ln -sf "$(realpath {input.R2})" {output.R2} >> {log} 2>&1
            """

################# 02 BEFORE QC #######################

rule before_fastqc:
    input:
        lambda wildcards: [
            f"{SCRATCH_DIR}/Concatenated/{wildcards.sample}_R1_concat.fastq.gz",
            f"{SCRATCH_DIR}/Concatenated/{wildcards.sample}_R2_concat.fastq.gz"
        ],
    output:
        html1=temp(f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_Before/{{sample}}_R1_concat_fastqc.html"),
        html2=temp(f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_Before/{{sample}}_R2_concat_fastqc.html"),
        zip1=f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_Before/{{sample}}_R1_concat_fastqc.zip",
        zip2=f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_Before/{{sample}}_R2_concat_fastqc.zip"
    conda:
        "../envs/QC_env.yaml"
    threads: config["fastqc"]["threads"]
    log:
        f"{LOG_DIR}/01_QC_Preprocessing/Before_QC/{{sample}}_QC.log"
    benchmark:
        f"{BENCH_DIR}/01_QC_Preprocessing/Before_QC/{{sample}}_QC.tsv"
    resources:
        mem_mb = config["fastqc"]["memory_mb"],
        runtime = config["fastqc"]["runtime_min"],
        cpus_per_task = config["fastqc"]["threads"]
    shell:
        """
        mkdir -p {OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_Before
        fastqc -t {threads} -o {OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_Before {input} > {log} 2>&1
        """

rule before_multiqc:
    input:
        expand(
            f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_Before/{{sample}}_R{{read}}_concat_fastqc.html",
            sample=SAMPLES,
            read=[1, 2]
        )
    output:
        outdir = directory(f"{OUTPUT_DIR}/02_Results/01_QC_Preprocessing/QC_Before")
    conda:
        "../envs/QC_env.yaml"
    threads: config["multiqc"]["threads"]
    log:
        f"{LOG_DIR}/01_QC_Preprocessing/Before_QC/Multi_QC.log"
    benchmark:
        f"{BENCH_DIR}/01_QC_Preprocessing/Before_QC/Multi_QC.tsv"
    resources:
        mem_mb = config["multiqc"]["memory_mb"],
        runtime = config["multiqc"]["runtime_min"],
        cpus_per_task = config["multiqc"]["threads"]
    shell:
        """
        mkdir -p {OUTPUT_DIR}/02_Results/01_QC_Preprocessing/QC_Before
        multiqc -o {output.outdir} {OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_Before > {log} 2>&1
        """

############### 03 TRIMMING #############################

rule Raw_Trimming:
    input:
        # ancient(): defense-in-depth so a merely-re-touched (not content-changed)
        # Concatenated/ file never spuriously triggers a rerun here — see
        # docs/log_Bug_fix.txt for why this link matters.
        R1=ancient(f"{SCRATCH_DIR}/Concatenated/{{sample}}_R1_concat.fastq.gz"),
        R2=ancient(f"{SCRATCH_DIR}/Concatenated/{{sample}}_R2_concat.fastq.gz"),
    output:
        # R1/R2_paired not temp() — see concat_lanes above for why (same cascade risk,
        # one step further down the chain, directly upstream of host_filtering).
        R1_paired=f"{SCRATCH_DIR}/Trimmed/{{sample}}_R1_paired.fastq.gz",
        R2_paired=f"{SCRATCH_DIR}/Trimmed/{{sample}}_R2_paired.fastq.gz",
        R1_unpaired=temp(f"{SCRATCH_DIR}/Trimmed/{{sample}}_R1_unpaired.fastq.gz"),
        R2_unpaired=temp(f"{SCRATCH_DIR}/Trimmed/{{sample}}_R2_unpaired.fastq.gz")
    threads: config["trimmomatic"]["threads"]
    params:
        nextera = config["trimmomatic"]["trimmomatic_adapters"],
        trailing= config["trimmomatic"]["trailing"],
        leading= config["trimmomatic"]["leading"],
        min_length= config["trimmomatic"]["minlen"]
    log:
        f"{LOG_DIR}/01_QC_Preprocessing/Trimming/{{sample}}_trimmomatic.log"
    benchmark:
        f"{BENCH_DIR}/01_QC_Preprocessing/Trimming/{{sample}}_trimmomatic.tsv"
    conda:
        "../envs/Trimmomatic.yaml"
    resources:
        mem_mb = config["trimmomatic"]["memory_mb"],
        runtime = config["trimmomatic"]["runtime_min"],
        cpus_per_task = config["trimmomatic"]["threads"]
    shell:
        "trimmomatic PE -phred33 -summary {log} -threads {threads} {input.R1} {input.R2} \
        {output.R1_paired} {output.R1_unpaired} {output.R2_paired} {output.R2_unpaired} \
        ILLUMINACLIP:{params.nextera}:2:30:10 \
        LEADING:{params.leading} TRAILING:{params.trailing} MINLEN:{params.min_length}"

############### 04 AFTER QC #####################################

rule After_fastqc:
    input:
        lambda wildcards: [
            f"{SCRATCH_DIR}/Trimmed/{wildcards.sample}_R1_paired.fastq.gz",
            f"{SCRATCH_DIR}/Trimmed/{wildcards.sample}_R2_paired.fastq.gz"
        ],
    output:
        html1=temp(f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_After/{{sample}}_R1_paired_fastqc.html"),
        html2=temp(f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_After/{{sample}}_R2_paired_fastqc.html"),
        zip1=f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_After/{{sample}}_R1_paired_fastqc.zip",
        zip2=f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_After/{{sample}}_R2_paired_fastqc.zip",
    conda:
        "../envs/QC_env.yaml"
    threads: config["fastqc"]["threads"]
    log:
        f"{LOG_DIR}/01_QC_Preprocessing/After_QC/{{sample}}_QC.log"
    benchmark:
        f"{BENCH_DIR}/01_QC_Preprocessing/After_QC/{{sample}}_QC.tsv"
    resources:
        mem_mb = config["fastqc"]["memory_mb"],
        runtime = config["fastqc"]["runtime_min"],
        cpus_per_task = config["fastqc"]["threads"]
    shell:
        """
        mkdir -p {OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_After
        fastqc -t {threads} -o {OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_After {input} > {log} 2>&1
        """


rule After_multiqc:
    input:
        expand(
            f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_After/{{sample}}_R{{read}}_paired_fastqc.html",
            sample=SAMPLES,
            read=[1, 2]
        )
    output:
        outdir = directory(f"{OUTPUT_DIR}/02_Results/01_QC_Preprocessing/QC_After")
    conda:
        "../envs/QC_env.yaml"
    threads: config["multiqc"]["threads"]
    log:
        f"{LOG_DIR}/01_QC_Preprocessing/After_QC/Multi_QC.log"
    benchmark:
        f"{BENCH_DIR}/01_QC_Preprocessing/After_QC/Multi_QC.tsv"
    resources:
        mem_mb = config["multiqc"]["memory_mb"],
        runtime = config["multiqc"]["runtime_min"],
        cpus_per_task = config["multiqc"]["threads"]
    shell:
        """
        mkdir -p {OUTPUT_DIR}/02_Results/01_QC_Preprocessing/QC_After
        multiqc -o {output.outdir} {OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/QC_After > {log} 2>&1
        """

###################### 05 TAX PROFILING ######################################

rule build_Kraken2_db:
    output:
        db = directory(config["kraken2_db"]["DB_dir"])
    params:
        host_list=lambda wildcards: " ".join(config["kraken2_db"]["refseq"])
    threads: config["kraken2_db"]["threads"]
    conda:
        "../envs/Kraken2.yaml"
    resources:
        mem_mb = config["kraken2_db"]["memory_mb"],
        runtime = config["kraken2_db"]["runtime_min"],
        disk_mb = config["kraken2_db"]["disk_mb"],
        cpus_per_task = config["kraken2_db"]["threads"]
    log:
        f"{LOG_DIR}/01_QC_Preprocessing/Kraken2/kraken2_db.log"
    benchmark:
        f"{BENCH_DIR}/01_QC_Preprocessing/Kraken2/kraken2_db.tsv"
    shell:
        """
        # Create the db folder
        mkdir -p $TMPDIR/kraken2_db
        # Build the database
        kraken2-build --download-taxonomy --db $TMPDIR/kraken2_db
        kraken2-build --download-library bacteria --db $TMPDIR/kraken2_db
        kraken2-build --download-library human --db $TMPDIR/kraken2_db
        echo "Adding files: ${params.host_list}" >> {log}
        for f in ${params.host_list}; do
            echo "File: $f" >> {log}
            kraken2-build --add-to-library "$f" --db $TMPDIR/kraken2_db
        done
        kraken2-build --build --threads {threads} --db $TMPDIR/kraken2_db

        # Copy the database back to the shared folder
        mkdir -p {output.db}
        cp -r $TMPDIR/kraken2_db/* {output.db}/
        """


rule run_kraken2:
    input:
        R1=f"{SCRATCH_DIR}/Trimmed/{{sample}}_R1_paired.fastq.gz",
        R2=f"{SCRATCH_DIR}/Trimmed/{{sample}}_R2_paired.fastq.gz",
        db= config["kraken2_db"]["DB_dir"]
    output:
        tab=temp(f"{SCRATCH_DIR}/Kraken2/{{sample}}_kraken2.out"),
        rep=temp(f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/Kraken2/{{sample}}_kraken2.report")
    conda:
        "../envs/Kraken2.yaml"
    threads: config["kraken2_run"]["threads"]
    log:
        f"{LOG_DIR}/01_QC_Preprocessing/Kraken2/{{sample}}_kraken2_run.log"
    benchmark:
        f"{BENCH_DIR}/01_QC_Preprocessing/Kraken2/{{sample}}_kraken2_run.tsv"
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


rule parse_kraken2_report:
    input:
        expand(f"{OUTPUT_DIR}/01_Analysis/01_QC_Preprocessing/Kraken2/{{sample}}_kraken2.report", sample=SAMPLES)
    output:
        f"{OUTPUT_DIR}/02_Results/01_QC_Preprocessing/Kraken2/Combined_kraken2_report.tsv"
    localrule: True
    conda:
        "../envs/Parsing.yaml"
    threads: 1
    log:
        f"{LOG_DIR}/01_QC_Preprocessing/Kraken2/Combined_kraken2_parsing.log"
    benchmark:
        f"{BENCH_DIR}/01_QC_Preprocessing/Kraken2/Combined_kraken2_parsing.tsv"
    shell:
        "python3 scripts/Combine_kraken_reports.py -o {output} -i {input} > {log} 2>&1"


########################## 06 HOST FILTERING ########################################

rule host_indexing:
    output:
        _HOST_INDEX_FILES
    conda:
        "../envs/bowtie2_mapping.yaml"
    threads: config["host_indexing"]["threads"]
    params:
        refs_seq    = lambda wildcards: " ".join(config["host_indexing"]["Host_refs"]),
        index_dir   = _HOST_INDEX_DIR,
        index_prefix = _HOST_INDEX_PREFIX,
    log:
        f"{LOG_DIR}/01_QC_Preprocessing/HF/host_indexing.log"
    benchmark:
        f"{BENCH_DIR}/01_QC_Preprocessing/HF/host_indexing.tsv"
    resources:
        mem_mb = config["host_indexing"]["memory_mb"],
        runtime = config["host_indexing"]["runtime_min"],
        cpus_per_task = config["host_indexing"]["threads"]
    shell:
        """
        cat {params.refs_seq} > $TMPDIR/hosts_combined.fa

        mkdir -p $TMPDIR/hosts_index
        bowtie2-build --threads {threads} $TMPDIR/hosts_combined.fa $TMPDIR/hosts_index/host_index > {log} 2>&1

        mkdir -p {params.index_dir}
        cp -r $TMPDIR/hosts_index/host_index.* {params.index_dir}/
        """

rule host_filtering:
    input:
        # ancient(): defense-in-depth against a spuriously fresher Trimmed/ mtime
        # cascading into a rerun of this rule and, transitively, assembly
        # (~23h/sample) — see docs/log_Bug_fix.txt.
        R1 = ancient(f"{SCRATCH_DIR}/Trimmed/{{sample}}_R1_paired.fastq.gz"),
        R2 = ancient(f"{SCRATCH_DIR}/Trimmed/{{sample}}_R2_paired.fastq.gz"),
        index = _HOST_INDEX_FILES
    output:
        R1_clean = f"{OUTPUT_DIR}/00_Data/02_Clean_Data/{{sample}}_R1_cleaned.fastq.gz",
        R2_clean = f"{OUTPUT_DIR}/00_Data/02_Clean_Data/{{sample}}_R2_cleaned.fastq.gz"
    conda:
        "../envs/bowtie2_mapping.yaml"
    threads: config["host_filtering"]["threads"]
    params:
        # Node-local scratch, not SCRATCH_DIR: this SAM file is never read again,
        # and temp() has no effect on params (only on declared outputs), so it
        # must not be left to accumulate on shared/persistent storage.
        sam = lambda wc: f"$TMPDIR/{wc.sample}_host.sam",
        Host_filt = f"{OUTPUT_DIR}/00_Data/02_Clean_Data/{{sample}}_R%_cleaned.fastq.gz",
        index = config["host_filtering"]["index_path"]
    log:
        f"{LOG_DIR}/01_QC_Preprocessing/HF/{{sample}}_host_filtering.log"
    benchmark:
        f"{BENCH_DIR}/01_QC_Preprocessing/HF/{{sample}}_host_filtering.tsv"
    resources:
        mem_mb = config["host_filtering"]["memory_mb"],
        runtime = config["host_filtering"]["runtime_min"],
        cpus_per_task = config["host_filtering"]["threads"]
    shell:
        """
        bowtie2 --very-sensitive-local --threads {threads} \
            --un-conc-gz {params.Host_filt} \
            -x {params.index} \
            -1 {input.R1} \
            -2 {input.R2} \
            -S {params.sam} 2> {log} \
        """
