# COMMON USED FUNCTION SNAKEFILE
# Last Updated: 2026-06-06
# Author: Matias Becker Burgos (Matias.BeckerBurgos@unil.ch)

#######################################################
# DESCRIPTION
#######################################################
# File containing commonly used functions in the pipelines
#
# Functions included:
#           - detect_samples()       : detect sample names from input FASTQ files
#           - get_excluded_samples() : compute blank/control exclusion list
#


##################################################
# IMPORT
##################################################

import os
import glob
import re
import sys
import pandas as pd
from snakemake.utils import *

##################################################
# FUNCTIONS
##################################################

def detect_samples(data_dir, pattern="*.fastq.gz"):
    """Detect unique sample names from FASTQ files in data_dir.

    Uses config['sample_pattern'] regex (group 1 = sample name).
    Silently skips files that don't match. Raises ValueError if no samples found.
    """
    files = glob.glob(os.path.join(data_dir, pattern))
    samples = set()
    regex = config.get("sample_pattern", r'^(.+?)_R[12].*\.fastq\.gz$')
    for f in files:
        base = os.path.basename(f)
        m = re.match(regex, base)
        if m is None:
            continue
        samples.add(m.group(1))
    if not samples:
        raise ValueError(
            f"detect_samples: no samples found in '{data_dir}' "
            f"matching pattern '{regex}'. "
            f"Check input_dir and sample_pattern in config.yaml."
        )
    return sorted(samples)


def get_excluded_samples(all_samples):
    """Return the union of explicitly listed and pattern-matched excluded samples.

    Reads config keys:
      exclude_samples : list of explicit sample names
      exclude_pattern : regex; empty string = disabled
    """
    excluded = set(config.get("exclude_samples", []))
    pattern = config.get("exclude_pattern", "")
    if pattern:
        excluded |= {s for s in all_samples if re.search(pattern, s)}
    return sorted(excluded)
