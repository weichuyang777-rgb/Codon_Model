#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
base_composition.py
===================
Teaching module, STEP 3 -- tally **base composition** of the human exon
coding region separately for each codon position.

Reproduces Table 1 of the PLOS Computational Biology Education manuscript
(counts in nucleotides):
    Position 1 : A 8487370  C 7648835  G 9762283  T 5338323
    Position 2 : A 9917738  C 7365391  G 5873697  T 8079985
    Position 3 : A 6418107  C 8877609  G 8700406  T 7240689

Standard library only.  The table is embedded as a constant so the script runs
standalone and reproduces the worked example.  To use real data, replace
BASE_COMP with counts parsed from a coding-sequence (CDS) file -- see README.md.

Run:
    python3 base_composition.py
"""

BASES = ['A', 'C', 'G', 'T']

# Embedded summary table: BASE_COMP[position][base] = nucleotide count.
# Values are taken directly from the manuscript's Table 1.
BASE_COMP = {
    1: {'A': 8487370, 'C': 7648835, 'G': 9762283, 'T': 5338323},
    2: {'A': 9917738, 'C': 7365391, 'G': 5873697, 'T': 8079985},
    3: {'A': 6418107, 'C': 8877609, 'G': 8700406, 'T': 7240689},
}


def main():
    print("=" * 64)
    print("Table 1. Base composition of the human exon coding region")
    print("=" * 64)
    header = "Position | " + " | ".join(f"{b:>10}" for b in BASES) + " |      Total"
    print(header)
    print("-" * len(header))
    grand = {b: 0 for b in BASES}
    for p in (1, 2, 3):
        row = f"   {p:<5}  | " + " | ".join(f"{BASE_COMP[p][b]:>10}" for b in BASES)
        total = sum(BASE_COMP[p].values())
        print(row + f" | {total:>10}")
        for b in BASES:
            grand[b] += BASE_COMP[p][b]
    total_line = f"   Total | " + " | ".join(f"{grand[b]:>10}" for b in BASES)
    print("-" * len(header))
    print(total_line + f" | {sum(grand.values()):>10}")

    # ---- Worked-results sentence: per-position richest / poorest base ------
    print("\nObservations (worked results):")
    for p in (1, 2, 3):
        richest = max(BASES, key=lambda b: BASE_COMP[p][b])
        poorest = min(BASES, key=lambda b: BASE_COMP[p][b])
        print(f"  Position {p}: richest base = {richest} "
              f"({BASE_COMP[p][richest]:,}); "
              f"poorest = {poorest} ({BASE_COMP[p][poorest]:,})")


if __name__ == '__main__':
    main()
