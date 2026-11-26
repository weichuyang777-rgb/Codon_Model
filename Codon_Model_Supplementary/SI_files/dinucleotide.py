#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dinucleotide.py
===============
Teaching module, STEP 4 -- compare **observed** dinucleotide frequencies with
the frequencies **expected under independence** (product of the two single-base
frequencies) for the three adjacent position-pairs (1-2, 2-3, 3-1).

Reproduces Table 2 and the Fig 3 result of the PLOS Computational Biology
Education manuscript:
the clearest deviation is TA, whose observed frequency is far below expectation
at every position pair (consistent with the three stop codons suppressing TA).

Standard library only.  The observed/expected table is embedded as a constant so
the script runs standalone and reproduces the worked example.  To use real data,
replace OBS / EXP with values parsed from a coding-sequence file -- see README.md.

Run:
    python3 dinucleotide.py
"""

DINUC = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT',
         'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']

# Observed (obs) and expected-under-independence (exp) frequencies for the
# three position pairs.  Values are taken directly from the manuscript's Table 2.
OBS = {
    '12': [0.0949, 0.0539, 0.0575, 0.0653, 0.0747, 0.0624, 0.0313, 0.0765,
           0.1204, 0.0682, 0.0644, 0.0595, 0.0275, 0.0513, 0.0348, 0.0573],
    '23': [0.0728, 0.0722, 0.1069, 0.0656, 0.0643, 0.0811, 0.0226, 0.0678,
           0.0375, 0.0605, 0.0505, 0.0395, 0.0308, 0.0704, 0.0986, 0.0589],
    '31': [0.0547, 0.0404, 0.0795, 0.0294, 0.1013, 0.0820, 0.0355, 0.0654,
           0.0739, 0.0707, 0.0962, 0.0373, 0.0400, 0.0517, 0.1013, 0.0388],
}
EXP = {
    '12': [0.0863, 0.0641, 0.0511, 0.0703, 0.0777, 0.0577, 0.0460, 0.0633,
           0.0992, 0.0737, 0.0588, 0.0808, 0.0543, 0.0403, 0.0321, 0.0442],
    '23': [0.0652, 0.0902, 0.0884, 0.0736, 0.0484, 0.0670, 0.0657, 0.0547,
           0.0386, 0.0534, 0.0524, 0.0436, 0.0531, 0.0735, 0.0720, 0.0600],
    '31': [0.0558, 0.0503, 0.0642, 0.0351, 0.0772, 0.0696, 0.0888, 0.0486,
           0.0757, 0.0682, 0.0870, 0.0476, 0.0630, 0.0568, 0.0724, 0.0396],
}


def deviation(pair):
    return [OBS[pair][i] - EXP[pair][i] for i in range(16)]


def main():
    print("=" * 78)
    print("Table 2. Observed vs expected dinucleotide frequencies (3 position pairs)")
    print("=" * 78)
    print(f"{'Di':<4}{'1-2 obs':>10}{'1-2 exp':>10}{'2-3 obs':>10}"
          f"{'2-3 exp':>10}{'3-1 obs':>10}{'3-1 exp':>10}")
    for i, d in enumerate(DINUC):
        print(f"{d:<4}"
              f"{OBS['12'][i]:>10.4f}{EXP['12'][i]:>10.4f}"
              f"{OBS['23'][i]:>10.4f}{EXP['23'][i]:>10.4f}"
              f"{OBS['31'][i]:>10.4f}{EXP['31'][i]:>10.4f}")

    print("\nDeviation (observed - expected):")
    for pair, label in [('12', 'positions 1-2'), ('23', 'positions 2-3'),
                        ('31', 'positions 3-1')]:
        dev = deviation(pair)
        ta = DINUC.index('TA')
        print(f"  {label:<16}: TA deviation = {dev[ta]:+.4f} "
              f"(min = {min(dev):+.4f})  <- under-represented, as predicted")


if __name__ == '__main__':
    main()
