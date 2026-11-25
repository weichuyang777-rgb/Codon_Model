#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
amino_acid_frequency.py
=======================
Teaching module, STEP 2 -- relate **synonymous-codon count** to **amino-acid
frequency in the human proteome** by ordinary least-squares linear regression.

Reproduces the regression in "Worked results and what they teach" of the
PLOS Computational Biology Education manuscript:
    Y(%) = 0.984 * X + 2.13 ,  R^2 = 0.51 ,  slope t = 4.33 ,  p < 0.001
and the grouped (by codon-count) fit, R^2 = 0.87.

Standard library only.  The 20 amino acids, their synonymous-codon count, and
their human-proteome frequency (%) are embedded as constants so the script runs
standalone and reproduces the worked example.  To use real data, replace AA/CC/FR
with values counted from a protein FASTA file (see README.md).

Run:
    python3 amino_acid_frequency.py
"""

from collections import defaultdict
import math


# Amino acids (single-letter), in the order used by the manuscript's Figure 1.
AA = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I',
      'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

# Number of synonymous codons encoding each amino acid (1..6).
CC = [4, 6, 2, 2, 2, 2, 2, 4, 2, 3, 6, 2, 1, 2, 4, 6, 4, 1, 2, 4]

# Observed frequency of each amino acid in the human proteome (percent).
FR = [8.4, 5.5, 4.1, 5.5, 2.2, 7.3, 4.2, 7.0, 2.3, 5.5,
      9.6, 5.7, 2.4, 4.0, 4.8, 6.9, 5.9, 1.2, 3.2, 6.9]


def ols(x, y):
    """Ordinary least squares. Returns (slope, intercept, r2, t)."""
    n = len(x)
    mx = sum(x) / n
    my = sum(y) / n
    sxx = sum((xi - mx) ** 2 for xi in x)
    sxy = sum((xi - mx) * (yi - my) for xi, yi in zip(x, y))
    syy = sum((yi - my) ** 2 for yi in y)
    slope = sxy / sxx
    intercept = my - slope * mx
    yhat = [slope * xi + intercept for xi in x]
    ss_res = sum((yi - yh) ** 2 for yi, yh in zip(y, yhat))
    r2 = 1 - ss_res / syy
    se_slope = math.sqrt(ss_res / ((n - 2) * sxx))
    t = slope / se_slope
    return slope, intercept, r2, t


def main():
    slope, intercept, r2, t = ols(CC, FR)
    print("=" * 64)
    print("Synonymous-codon count  vs  proteome frequency")
    print("=" * 64)
    print(f"  OLS fit :  Y(%) = {slope:.3f} * X + {intercept:.2f}")
    print(f"  R^2     = {r2:.2f}")
    print(f"  slope t = {t:.2f}   (p < 0.001)")
    print(f"  proportion form:  Y = {slope/100:.5f} * X + {intercept/100:.4f}")

    # Group amino acids by codon count, average the frequency, refit.
    groups = defaultdict(list)
    for x, y in zip(CC, FR):
        groups[x].append(y)
    gx = sorted(groups)
    gy = [sum(groups[k]) / len(groups[k]) for k in gx]
    _, _, r2g, _ = ols(gx, gy)
    print("\n  Grouped by codon-count (5 group means):")
    print(f"    codon-count groups : {gx}")
    print(f"    mean frequencies   : {[round(v, 2) for v in gy]}")
    print(f"    R^2 = {r2g:.2f}   <- teaching point: aggregation inflates R^2")

    print("\nInterpretation: more codons -> higher average frequency, but the")
    print("moderate R^2 shows the relationship is far from deterministic;" )
    print("residual variation is driven by codon usage bias.")


if __name__ == '__main__':
    main()
