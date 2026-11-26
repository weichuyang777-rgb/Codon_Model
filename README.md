# Codon_Model

This repository contains the manuscript, cover letter, and supplementary reproducibility package for the *PLOS Computational Biology* (Education) submission:

> **"A Python-based computational module for teaching codon degeneracy and codon usage bias"**  
> Yang Lei — School of Laboratory Medicine and Life Sciences, Wenzhou Medical University, Wenzhou, Zhejiang, China

## Purpose

The project turns the standard genetic code into a hands-on computational biology tutorial consisting of four Python exercises:

1. Encoding the genetic code as a nested Python dictionary and counting positional degeneracy (`codon_model.py`).
2. Testing positional base differences using paired *t*-tests (`codon_model.py`).
3. Regressing human proteome amino-acid frequencies on synonymous codon counts (`amino_acid_frequency.py`).
4. Analyzing per-position base composition and dinucleotide frequencies in human coding exons (`base_composition.py`, `dinucleotide.py`).

All analysis and figure-generation scripts use **only the Python 3 standard library** (no third-party packages required) and reproduce every number, table, and SVG figure reported in the manuscript.

## Repository Structure

```text
Codon_Model/
├── Codon_Model_PLOS_Comput_Biol_Education.docx   # Main manuscript (PLOS Computational Biology Education)
├── Cover_Letter_PLOS.docx                        # Submission cover letter (.docx)
├── Cover_Letter_PLOS.txt                         # Submission cover letter (plain text)
└── Codon_Model_Supplementary/                    # Supplementary code, figures, and metadata
    ├── codon_model.py                            # Genetic-code dictionary, degeneracy counts, paired t-tests
    ├── amino_acid_frequency.py                   # OLS regression of codon count vs. proteome frequency
    ├── base_composition.py                       # Per-position base composition (Table 1)
    ├── dinucleotide.py                           # Observed vs. expected dinucleotide frequencies (Table 2)
    ├── make_figures.py                           # Generates Figs. 1–3 as SVG (pure stdlib)
    ├── Fig1_codon_frequency.svg                  # Figure 1
    ├── Fig2_base_composition.svg                 # Figure 2
    ├── Fig3_dinucleotide_deviation.svg           # Figure 3
    ├── CITATION.cff                              # Citation metadata
    ├── LICENSE                                   # MIT License
    └── README.md                                 # Detailed supplementary documentation
```

## Quick Start

```bash
cd Codon_Model_Supplementary
python3 codon_model.py
python3 amino_acid_frequency.py
python3 base_composition.py
python3 dinucleotide.py
python3 make_figures.py
```

For detailed expected outputs and instructions on re-running the pipeline with raw UCSC Genome Browser tables, see [Codon_Model_Supplementary/README.md](Codon_Model_Supplementary/README.md).

## Citation & License

- **License**: MIT License (see [Codon_Model_Supplementary/LICENSE](Codon_Model_Supplementary/LICENSE)).
- **Citation**: See [Codon_Model_Supplementary/CITATION.cff](Codon_Model_Supplementary/CITATION.cff).
