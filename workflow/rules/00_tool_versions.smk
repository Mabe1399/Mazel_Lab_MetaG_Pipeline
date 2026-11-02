# TOOL VERSION PROBES SNAKEFILE
# Last Updated: 2026-11-02
# Author: Matias Becker Burgos (Matias.BeckerBurgos@unil.ch)

#######################################################
# DESCRIPTION
#######################################################
# MODULE INTERFACE
#   INPUT:     —
#   OUTPUT:    {OUTPUT_DIR}/01_Analysis/ToolVersions/{env}_versions.tsv  (tool<TAB>version, one file per env)
#   REQUIRES:  —
#   STEPS KEY: — (always included; _get_tool_version_inputs() in the Snakefile requests only the
#              probes of envs that the enabled steps use, and aggregate_tool_versions merges them)
#   USES:      —
#   NOTE:      One probe per conda env, kept in one file because the grouped envs span several
#              modules (Core: 01/01b/02/02b + every Python parsing rule + 05; Binning: 03/03b).
#              Adding a tool to an env means adding a line to its probe here. Same layout as
#              03_Community_Analysis_Pipeline's rules/00_tool_versions.smk.

# Version strings use the same extraction everywhere: several tools print banners/warnings
# before the version line, so a naive first-line parse isn't reliable. -a: bowtie2's version
# banner contains non-text bytes, and without it grep prints "Binary file matches" instead.
_VERSION_RE = "grep -aoE '[0-9]+\\.[0-9]+(\\.[0-9]+)?' | head -1 || echo unknown"


rule tool_version_core:
    output:
        versions = temp(f"{OUTPUT_DIR}/01_Analysis/ToolVersions/core_versions.tsv")
    params:
        v = _VERSION_RE
    conda:
        "../envs/Core.yaml"
    localrule: True
    shell:
        """
        mkdir -p $(dirname {output.versions})
        {{
            printf "fastqc\\t%s\\n" "$(fastqc --version 2>&1 | {params.v})"
            printf "multiqc\\t%s\\n" "$(multiqc --version 2>&1 | {params.v})"
            printf "trimmomatic\\t%s\\n" "$(trimmomatic -version 2>&1 | {params.v})"
            printf "bowtie2\\t%s\\n" "$(bowtie2 --version 2>&1 | {params.v})"
            printf "kraken2\\t%s\\n" "$(kraken2 --version 2>&1 | {params.v})"
            printf "spades\\t%s\\n" "$(spades.py --version 2>&1 | {params.v})"
            printf "quast\\t%s\\n" "$(quast.py --version 2>&1 | {params.v})"
            printf "pandas\\t%s\\n" "$(python3 -c 'import pandas; print(pandas.__version__)' 2>&1 | {params.v})"
        }} > {output.versions}
        """


rule tool_version_binning:
    output:
        versions = temp(f"{OUTPUT_DIR}/01_Analysis/ToolVersions/binning_versions.tsv")
    params:
        v = _VERSION_RE
    conda:
        "../envs/Binning.yaml"
    localrule: True
    shell:
        """
        mkdir -p $(dirname {output.versions})
        {{
            printf "fairy\\t%s\\n" "$(fairy --version 2>&1 | {params.v})"
            printf "metabat2\\t%s\\n" "$(metabat2 --help 2>&1 | {params.v})"
            printf "maxbin2\\t%s\\n" "$(run_MaxBin.pl -v 2>&1 | {params.v})"
            printf "das_tool\\t%s\\n" "$(DAS_Tool --version 2>&1 | {params.v})"
            printf "checkm2\\t%s\\n" "$(checkm2 --version 2>&1 | {params.v})"
        }} > {output.versions}
        """


rule tool_version_semibin2:
    output:
        versions = temp(f"{OUTPUT_DIR}/01_Analysis/ToolVersions/semibin2_versions.tsv")
    params:
        v = _VERSION_RE
    conda:
        "../envs/SemiBin2.yaml"
    localrule: True
    shell:
        """
        mkdir -p $(dirname {output.versions})
        printf "semibin2\\t%s\\n" "$(SemiBin2 --version 2>&1 | {params.v})" > {output.versions}
        """


rule tool_version_reporting:
    output:
        versions = temp(f"{OUTPUT_DIR}/01_Analysis/ToolVersions/reporting_versions.tsv")
    params:
        v = _VERSION_RE
    conda:
        "../envs/Reporting.yaml"
    localrule: True
    shell:
        """
        mkdir -p $(dirname {output.versions})
        {{
            printf "gtdbtk\\t%s\\n" "$(gtdbtk --version 2>&1 | {params.v})"
            printf "drep\\t%s\\n" "$(dRep -h 2>&1 | {params.v})"
        }} > {output.versions}
        """
