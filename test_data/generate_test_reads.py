#!/usr/bin/env python3
"""Generate small synthetic paired-end FASTQ test data for cluster smoke-testing.

Simulates a 2-sample, 3-"species" mock community with differential relative
abundance between samples (so binning tools have a real coverage signal to
work with) and writes gzipped, pre-cleaned FASTQ pairs directly to reads/ —
matching the {sample}_R{1,2}_cleaned.fastq.gz contract expected when
steps.QC: false in config.yaml (see test_data/config.test.yaml).

Deterministic (fixed seed) — rerun this script to regenerate reads/ from
scratch rather than hand-editing the FASTQ files.

Usage: python3 generate_test_reads.py
"""

import gzip
import random

random.seed(42)

READ_LENGTH = 150
FRAGMENT_MEAN = 350
FRAGMENT_SD = 40
PER_BASE_ERROR_RATE = 0.005
BASES = "ACGT"

GENOMES = {
    "speciesA": 80_000,
    "speciesB": 60_000,
    "speciesC": 40_000,
}

# (sample_name, {species: coverage_depth})
SAMPLES = [
    ("sampleA", {"speciesA": 40, "speciesB": 15, "speciesC": 3}),
    ("sampleB", {"speciesA": 3, "speciesB": 40, "speciesC": 20}),
]


def random_genome(length):
    return "".join(random.choice(BASES) for _ in range(length))


def revcomp(seq):
    comp = str.maketrans("ACGT", "TGCA")
    return seq.translate(comp)[::-1]


def add_errors(seq, rate):
    seq = list(seq)
    for i in range(len(seq)):
        if random.random() < rate:
            seq[i] = random.choice([b for b in BASES if b != seq[i]])
    return "".join(seq)


def simulate_reads_for_genome(genome, depth, read_id_prefix):
    n_fragments = max(1, int(len(genome) * depth / (2 * READ_LENGTH)))
    records = []
    for i in range(n_fragments):
        frag_len = max(
            READ_LENGTH * 2, int(random.gauss(FRAGMENT_MEAN, FRAGMENT_SD))
        )
        if frag_len >= len(genome):
            continue
        start = random.randint(0, len(genome) - frag_len)
        fragment = genome[start : start + frag_len]
        r1 = add_errors(fragment[:READ_LENGTH], PER_BASE_ERROR_RATE)
        r2 = add_errors(revcomp(fragment[-READ_LENGTH:]), PER_BASE_ERROR_RATE)
        read_name = f"{read_id_prefix}_{i}"
        records.append((read_name, r1, r2))
    return records


def write_fastq_pair(records, r1_path, r2_path):
    with gzip.open(r1_path, "wt") as f1, gzip.open(r2_path, "wt") as f2:
        for name, r1, r2 in records:
            qual = "I" * READ_LENGTH
            f1.write(f"@{name}/1\n{r1}\n+\n{qual}\n")
            f2.write(f"@{name}/2\n{r2}\n+\n{qual}\n")


def main():
    genomes = {name: random_genome(length) for name, length in GENOMES.items()}

    for sample_name, abundances in SAMPLES:
        all_records = []
        for species, depth in abundances.items():
            all_records.extend(
                simulate_reads_for_genome(
                    genomes[species], depth, f"{sample_name}_{species}"
                )
            )
        random.shuffle(all_records)

        r1_path = f"reads/{sample_name}_R1_cleaned.fastq.gz"
        r2_path = f"reads/{sample_name}_R2_cleaned.fastq.gz"
        write_fastq_pair(all_records, r1_path, r2_path)
        print(f"{sample_name}: {len(all_records)} read pairs -> {r1_path}, {r2_path}")


if __name__ == "__main__":
    main()
