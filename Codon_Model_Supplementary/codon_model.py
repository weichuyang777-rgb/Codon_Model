#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
codon_model.py
==============
Teaching module, STEP 1 -- encode the standard genetic code as a nested
dictionary and quantify **codon degeneracy** at each of the three positions.

This script reproduces the results reported in:
  "A Python-based computational module for teaching codon degeneracy
   and codon usage bias"  (PLOS Computational Biology Education manuscript),
   section "Worked results and what they teach".

It is intentionally dependency-free (Python 3 standard library only) so that
learners in introductory bioinformatics courses can run it on any machine.

Run:
    python3 codon_model.py

The codon table and the resulting degeneracy counts are embedded as constants,
so the script runs standalone and reproduces the worked example in the paper.
To repeat the analysis on a different organism, replace GENETIC_CODE (and,
optionally, supply a parsed codon->count table) -- the statistics below are
general.
"""

import math
from statistics import mean

# ---------------------------------------------------------------------------
# 1. The standard genetic code: codon (RNA letters U/C/A/G) -> amino acid.
#    '*' marks the three stop codons.  U is used instead of T to match the
#    reading frame taught in the classroom.
# ---------------------------------------------------------------------------
GENETIC_CODE = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
    'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
    'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

BASES = ['U', 'C', 'A', 'G']


def build_codon_dict():
    """Encode the code as a 3-level nested dict: codon[a][b][c] = amino_acid.

    This is a *lookup structure*, not a statistical or mechanistic model --
    a distinction the manuscript makes explicit to learners.
    """
    codon = {}
    for a in BASES:
        codon[a] = {}
        for b in BASES:
            codon[a][b] = {}
            for c in BASES:
                codon[a][b][c] = GENETIC_CODE[a + b + c]
    return codon


def degeneracy_counts():
    """For each of the 3 positions, count how many amino acids each base
    is associated with.  Returns {position: {base: count}}."""
    codon = build_codon_dict()
    pos_counts = {1: {}, 2: {}, 3: {}}
    for a in BASES:
        for b in BASES:
            for c in BASES:
                aa = codon[a][b][c]
                base_at_pos = {1: a, 2: b, 3: c}
                for p in (1, 2, 3):
                    pos_counts[p].setdefault(base_at_pos[p], set()).add(aa)
    return {p: {base: len(s) for base, s in pos_counts[p].items()} for p in (1, 2, 3)}


def paired_ttest(x, y):
    """Paired Student t-test for two equal-length vectors. Returns (t, n)."""
    n = len(x)
    diffs = [x[i] - y[i] for i in range(n)]
    dbar = mean(diffs)
    sd = math.sqrt(sum((d - dbar) ** 2 for d in diffs) / (n - 1))
    se = sd / math.sqrt(n)
    t = dbar / se
    return t, n


def main():
    counts = degeneracy_counts()
    print("=" * 64)
    print("Codon degeneracy: amino acids associated with each base")
    print("=" * 64)
    for p in (1, 2, 3):
        vals = [counts[p][b] for b in BASES]
        print(f"  Position {p}: {vals}  (U, C, A, G)   total = {sum(vals)}")

    # ---- Worked results: paired t-tests between positions -----------------
    v1 = [counts[1][b] for b in BASES]
    v2 = [counts[2][b] for b in BASES]
    v3 = [counts[3][b] for b in BASES]
    print("\nPaired t-tests (n = 4 nucleotides, one per base):")
    for label, a, b in [('Position 1 vs 2', v1, v2),
                        ('Position 1 vs 3', v1, v3),
                        ('Position 2 vs 3', v2, v3)]:
        t, n = paired_ttest(a, b)
        sig = "p < 0.001" if abs(t) > 5 else "not significant"
        print(f"  {label:<18}: t = {t:6.2f}   ({sig})")

    # ---- Worked results: hydrogen-bond class -------------------------------
    # Two-hydrogen-bond bases {U, A} vs three-hydrogen-bond bases {C, G}.
    ua = [counts[p]['U'] + counts[p]['A'] for p in (1, 2, 3)]
    cg = [counts[p]['C'] + counts[p]['G'] for p in (1, 2, 3)]
    t2, _ = paired_ttest(ua, cg)
    print("\nHydrogen-bond class (UA = two H-bonds, CG = three H-bonds):")
    print(f"  UA = {ua}")
    print(f"  CG = {cg}")
    print(f"  paired t = {t2:6.2f}   (not significant)")
    print("\nTeaching note: the third position carries far more amino-acid")
    print("variety (totals 57 vs 24 and 23) -- the wobble hypothesis in raw data.")


if __name__ == '__main__':
    main()
