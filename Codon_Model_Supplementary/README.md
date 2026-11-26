# Supplementary Material — Codon Degeneracy & Codon Usage Bias Teaching Module

This folder contains the **Python scripts** that reproduce every number, table,
and figure in the teaching manuscript:

> *"A Python-based computational module for teaching codon degeneracy and
>  codon usage bias"* (submitted as an **Education** article to *PLOS Computational Biology*).

All scripts use **only the Python 3 standard library** — no numpy, no
matplotlib, no machine-learning package. They run on any machine with Python 3
installed and reproduce the worked example from the paper out of the box (the
key summary tables are embedded as constants). This keeps the manuscript's
reproducibility claim exact.

---

## File manifest

| File | Reproduces | Manuscript section |
|------|------------|--------------------|
| `codon_model.py` | Genetic-code dictionary, degeneracy counts, paired *t*-tests | Design of the module; Worked results |
| `amino_acid_frequency.py` | OLS regression of codon count vs proteome frequency | Worked results, Fig. 1 |
| `base_composition.py` | Per-position base composition table (Table 1) | Worked results, Table 1, Fig. 2 |
| `dinucleotide.py` | Observed vs expected dinucleotide frequencies (Table 2) | Worked results, Table 2, Fig. 3 |
| `make_figures.py` | The three figures as SVG (pure stdlib) | Figs. 1–3 |
| `Fig1_*.svg`, `Fig2_*.svg`, `Fig3_*.svg` | Generated figures | Figs. 1–3 |
| `CITATION.cff` | Machine-readable citation metadata | — |
| `LICENSE` | MIT license text | — |
| `README.md` | This file | — |

---

## How to run

```bash
python3 codon_model.py
python3 amino_acid_frequency.py
python3 base_composition.py
python3 dinucleotide.py
python3 make_figures.py        # writes Fig1/2/3_*.svg
```

Expected headline outputs (for quick verification):

- Degeneracy totals: Position 1 = 24, Position 2 = 23, Position 3 = 57.
- Paired *t*: Pos1–Pos2 = 0.33 (n.s.); Pos1–Pos3 = 9.66; Pos2–Pos3 = 6.43 (both *p* < 0.001).
- Regression: *Y*(%) = 0.984·*X* + 2.13, *R*² = 0.51, slope *t* = 4.33, *p* < 0.001.
- Grouped (by codon count) fit: *R*² = 0.87.
- Dinucleotide TA is the strongest under-representation at all three position pairs.

---

## Using real data (optional)

The embedded constants are the exact values reported in the paper. To repeat the
analysis on the raw UCSC downloads (or another organism):

1. Download from the **UCSC Genome Browser → Table Browser**
   (<http://genome.ucsc.edu/cgi-bin/hgTables>):
   - Human protein sequences (`hgTables_protein.txt`, ~28.6 MB).
   - Human exon / coding-sequence regions (`hgTables_exon.txt`, ~97.6 MB).
2. In `amino_acid_frequency.py`, replace `FR` with residue counts parsed from the
   protein FASTA.
3. In `base_composition.py` / `dinucleotide.py`, replace `BASE_COMP` / `OBS`+`EXP`
   with counts parsed from the coding-sequence file. The statistics and plotting
   code above remain unchanged.

No private or restricted data are used.

---

## Deposition and citation

The exact version of the code accompanying the manuscript is archived in the
public GitHub repository at
<https://github.com/weichuyang777-rgb/Codon_model/commit/4829f553eba8a598ddd2a9a1d5999db108a1a14e>.
This package reproduces that version. See `CITATION.cff` in this folder
for the machine-readable citation metadata.

No DOI is assigned for this deposit. The code is cited by the repository URL
together with the commit identifier
`4829f553eba8a598ddd2a9a1d5999db108a1a14e`, which pins the exact version used in
the manuscript; GitHub repositories can change over time, so the commit
identifier should always be included when retrieving or citing the code. A
zipped copy of this folder (`Codon_Model_Supplementary.zip`) is also provided
alongside the manuscript.

---

## License

This software is released under the **MIT License** (see `LICENSE`). It may be
freely used, modified and redistributed, including for teaching and commercial
purposes, provided that the copyright notice and this permission notice are
retained. The same license (`MIT`) is declared in `CITATION.cff`.
