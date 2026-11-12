#!/usr/bin/env python3
"""Build a small test contig set from a real assembly.

Takes a filtered metaSPAdes assembly (the Contig_filter output,
{sample}_min1500_contigs.fasta) and writes a deterministic, stratified subset:

- contigs are drawn per length class, in roughly the assembly's proportions,
  so the share of 1000-2500 bp sequence (which sets SemiBin2's minimum contig
  length: 2500 bp if it is >= 5%, else 1000 bp) stays the same as in the source;
- edge cases are always included: the shortest contig, contigs of exactly
  2500 bp (SemiBin2's cutoff) and 2501 bp, odd lengths (unequal halves after
  SemiBin2 split_contigs), the highest-coverage contig and the GC extremes;
- headers and sequences are copied unchanged, in the source's NODE order.

The source assembly is too large to commit (hundreds of MB); keep it next to
this script under contigs/ (gitignored) and rerun this script to regenerate the
subset rather than editing it by hand.

Usage:
    python3 make_test_contigs.py contigs/Glomus_China25_F_10_min1500_contigs.fasta \
        contigs/F10_subset_min1500_contigs.fasta
"""

import random
import re
import sys

SEED = 42

# (min length, max length (exclusive), number of contigs to draw)
LENGTH_CLASSES = [
    (1_500, 2_500, 600),
    (2_500, 5_000, 300),
    (5_000, 10_000, 120),
    (10_000, 50_000, 60),
    (50_000, 10**9, 6),
]

HEADER = re.compile(r"^NODE_(\d+)_length_(\d+)_cov_([\d.]+)$")
LINE_WIDTH = 60


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:], []
            else:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)


def gc(seq):
    return (seq.count("G") + seq.count("C")) / len(seq)


def short_fraction(lengths):
    """Share of bp in 1000-2500 bp contigs, as SemiBin2's load_fasta computes it."""
    return sum(n for n in lengths if 1000 <= n <= 2500) / sum(lengths)


def main(src, dst):
    contigs = []
    for name, seq in read_fasta(src):
        m = HEADER.match(name)
        if not m:
            sys.exit(f"Unexpected header (not metaSPAdes NODE_..._cov_...): {name}")
        if int(m.group(2)) != len(seq):
            sys.exit(f"Header length does not match sequence length: {name}")
        contigs.append({"name": name, "seq": seq, "node": int(m.group(1)),
                        "len": len(seq), "cov": float(m.group(3)), "gc": gc(seq)})

    # Edge cases, always kept
    picked = set()
    by_len = lambda n: [c["node"] for c in contigs if c["len"] == n]
    picked.add(min(contigs, key=lambda c: c["len"])["node"])
    picked.update(by_len(2500)[:2])
    picked.update(by_len(2501)[:2])
    picked.add(max(contigs, key=lambda c: c["cov"])["node"])
    picked.add(min(contigs, key=lambda c: c["gc"])["node"])
    picked.add(max(contigs, key=lambda c: c["gc"])["node"])

    # Stratified random draw per length class
    rng = random.Random(SEED)
    for lo, hi, n in LENGTH_CLASSES:
        pool = [c["node"] for c in contigs if lo <= c["len"] < hi and c["node"] not in picked]
        picked.update(rng.sample(pool, min(n, len(pool))))

    subset = [c for c in contigs if c["node"] in picked]  # source (NODE) order
    with open(dst, "w") as out:
        for c in subset:
            out.write(f">{c['name']}\n")
            for i in range(0, c["len"], LINE_WIDTH):
                out.write(c["seq"][i:i + LINE_WIDTH] + "\n")

    # Checks: the subset must behave like the source where the pipeline cares
    src_len = [c["len"] for c in contigs]
    sub_len = [c["len"] for c in subset]
    src_short, sub_short = short_fraction(src_len), short_fraction(sub_len)
    assert (src_short >= 0.05) == (sub_short >= 0.05), "SemiBin2 minimum length would differ from the source"
    assert len({c["name"] for c in subset}) == len(subset), "duplicate contig names"
    assert any(n % 2 for n in sub_len), "no odd-length contigs"

    print(f"source: {len(contigs):,} contigs, {sum(src_len) / 1e6:.1f} Mb, 1000-2500 bp share {src_short:.1%}")
    print(f"subset: {len(subset):,} contigs, {sum(sub_len) / 1e6:.2f} Mb, 1000-2500 bp share {sub_short:.1%}")
    print(f"        SemiBin2 minimum length: {2500 if sub_short >= 0.05 else 1000} bp, "
          f"{sum(1 for n in sub_len if n >= 2500):,} contigs binned by SemiBin2")
    print(f"        lengths {min(sub_len):,}-{max(sub_len):,} bp, "
          f"{sum(n % 2 for n in sub_len)} odd, {sum(1 for n in sub_len if n == 2500)} at exactly 2500 bp")
    print(f"wrote {dst}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
