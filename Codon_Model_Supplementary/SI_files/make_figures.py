#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_figures.py
===============
Generate the three manuscript figures as **SVG** using only the Python 3
standard library (no matplotlib / numpy required).  This keeps the manuscript's
reproducibility claim -- "Python 3 and a few tens of megabytes of public data,
no external package" -- exact.

If `matplotlib` is installed it is NOT used here on purpose: SVG is text, fully
reproducible, version-controllable and suitable for journal submission.  The
embedded data
arrays are identical to those used to draw the submitted PNG figures.

Outputs (in the current working directory):
    Fig1_codon_frequency.svg
    Fig2_base_composition.svg
    Fig3_dinucleotide_deviation.svg

Run:
    python3 make_figures.py
"""

# ----------------------------- Figure 1 data ------------------------------
AA = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I',
      'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
CC = [4, 6, 2, 2, 2, 2, 2, 4, 2, 3, 6, 2, 1, 2, 4, 6, 4, 1, 2, 4]
FR = [8.4, 5.5, 4.1, 5.5, 2.2, 7.3, 4.2, 7.0, 2.3, 5.5,
      9.6, 5.7, 2.4, 4.0, 4.8, 6.9, 5.9, 1.2, 3.2, 6.9]
SLOPE, INTER = 0.984, 2.13

# colour by synonymous-codon count
CC_COLORS = {1: '#800080', 2: '#FF0000', 3: '#008000', 4: '#0000FF', 6: '#FF8C00'}

# ----------------------------- Figure 2 data ------------------------------
POS = ['Pos 1', 'Pos 2', 'Pos 3']
BASES = ['A', 'C', 'G', 'T']
BASE_COLORS = {'A': '#228B22', 'C': '#4169E1', 'G': '#B8860B', 'T': '#DC143C'}
FIG2 = {
    'Pos 1': [8487370, 7648835, 9762283, 5338323],
    'Pos 2': [9917738, 7365391, 5873697, 8079985],
    'Pos 3': [6418107, 8877609, 8700406, 7240689],
}

# ----------------------------- Figure 3 data ------------------------------
DINUC = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT',
         'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
O12 = [0.0949, 0.0539, 0.0575, 0.0653, 0.0747, 0.0624, 0.0313, 0.0765,
       0.1204, 0.0682, 0.0644, 0.0595, 0.0275, 0.0513, 0.0348, 0.0573]
E12 = [0.0863, 0.0641, 0.0511, 0.0703, 0.0777, 0.0577, 0.0460, 0.0633,
       0.0992, 0.0737, 0.0588, 0.0808, 0.0543, 0.0403, 0.0321, 0.0442]
O23 = [0.0728, 0.0722, 0.1069, 0.0656, 0.0643, 0.0811, 0.0226, 0.0678,
       0.0375, 0.0605, 0.0505, 0.0395, 0.0308, 0.0704, 0.0986, 0.0589]
E23 = [0.0652, 0.0902, 0.0884, 0.0736, 0.0484, 0.0670, 0.0657, 0.0547,
       0.0386, 0.0534, 0.0524, 0.0436, 0.0531, 0.0735, 0.0720, 0.0600]
O31 = [0.0547, 0.0404, 0.0795, 0.0294, 0.1013, 0.0820, 0.0355, 0.0654,
       0.0739, 0.0707, 0.0962, 0.0373, 0.0400, 0.0517, 0.1013, 0.0388]
E31 = [0.0558, 0.0503, 0.0642, 0.0351, 0.0772, 0.0696, 0.0888, 0.0486,
       0.0757, 0.0682, 0.0870, 0.0476, 0.0630, 0.0568, 0.0724, 0.0396]


def svg_head(w, h, title):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" '
            f'viewBox="0 0 {w} {h}" font-family="Arial, sans-serif">'
            f'<rect width="{w}" height="{h}" fill="white"/>'
            f'<text x="50%" y="26" text-anchor="middle" font-size="15" '
            f'font-weight="bold">{title}</text>')


def fig1():
    w, h = 620, 440
    L, R, T, B = 75, 575, 70, 370
    s = [svg_head(w, h, "Fig. 1.  Synonymous-codon count versus proteome frequency")]
    # grid + y axis labels
    for yv in range(0, 12, 2):
        py = B - (yv / 11) * (B - T)
        s.append(f'<line x1="{L}" y1="{py}" x2="{R}" y2="{py}" stroke="#DDDDDD"/>')
        s.append(f'<text x="{L-8}" y="{py+4}" text-anchor="end" font-size="10">{yv}</text>')
    for xv in range(1, 7):
        px = L + ((xv - 1) / 5) * (R - L)
        s.append(f'<text x="{px}" y="{B+20}" text-anchor="middle" font-size="10">{xv}</text>')
    s.append(f'<text x="14" y="44" font-size="12">Frequency (%)</text>')
    s.append(f'<text x="{(L+R)/2}" y="{B+40}" text-anchor="middle" font-size="12">'
             f'Number of synonymous codons</text>')
    s.append(f'<line x1="{L}" y1="{T}" x2="{L}" y2="{B}" stroke="black"/>')
    s.append(f'<line x1="{L}" y1="{B}" x2="{R}" y2="{B}" stroke="black"/>')
    for i in range(20):
        px = L + ((CC[i] - 1) / 5) * (R - L)
        py = B - (FR[i] / 11) * (B - T)
        col = CC_COLORS[CC[i]]
        s.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="5" fill="{col}" '
                 f'stroke="black" stroke-width="1"/>')
    # regression line
    x1, x2 = 1, 6
    y1 = SLOPE * x1 + INTER
    y2 = SLOPE * x2 + INTER
    px1, px2 = L, R
    py1 = B - (y1 / 11) * (B - T)
    py2 = B - (y2 / 11) * (B - T)
    s.append(f'<line x1="{px1}" y1="{py1:.1f}" x2="{px2}" y2="{py2:.1f}" '
             f'stroke="black" stroke-width="2"/>')
    s.append(f'<text x="{L+12}" y="{T+14}" font-size="10">OLS fit: y = 0.984x + 2.13</text>')
    s.append(f'<text x="{L+12}" y="{T+30}" font-size="10">R2 = 0.51, p &lt; 0.001</text>')
    s.append('</svg>')
    return "\n".join(s)


def fig2():
    w, h = 640, 440
    L, R, T, B = 70, 600, 70, 370
    s = [svg_head(w, h, "Fig. 2.  Base composition across the three codon positions (human exons)")]
    for yv in range(0, 12, 2):
        py = B - (yv / 11) * (B - T)
        s.append(f'<line x1="{L}" y1="{py}" x2="{R}" y2="{py}" stroke="#DDDDDD"/>')
        s.append(f'<text x="{L-8}" y="{py+4}" text-anchor="end" font-size="10">{yv}</text>')
    s.append(f'<text x="8" y="44" font-size="12">Base count (millions)</text>')
    s.append(f'<text x="{(L+R)/2}" y="{B+40}" text-anchor="middle" font-size="12">Codon position</text>')
    s.append(f'<line x1="{L}" y1="{T}" x2="{L}" y2="{B}" stroke="black"/>')
    s.append(f'<line x1="{L}" y1="{B}" x2="{R}" y2="{B}" stroke="black"/>')
    gw = (R - L) / 3
    for pi in range(3):
        gx0 = L + pi * gw
        s.append(f'<text x="{gx0+gw/2}" y="{B+18}" text-anchor="middle" font-size="10">{POS[pi]}</text>')
        bw = (gw - 20) / 4
        for bi in range(4):
            val = FIG2[POS[pi]][bi] / 1e6
            bh = (val / 11) * (B - T)
            bx = gx0 + 10 + bi * bw
            by = B - bh
            s.append(f'<rect x="{bx:.1f}" y="{by:.1f}" width="{bw-3:.1f}" height="{bh:.1f}" '
                     f'fill="{BASE_COLORS[BASES[bi]]}"/>')
    # legend
    lx, ly = R - 120, T + 6
    for bi in range(4):
        s.append(f'<rect x="{lx}" y="{ly+bi*16}" width="12" height="12" fill="{BASE_COLORS[BASES[bi]]}"/>')
        s.append(f'<text x="{lx+18}" y="{ly+bi*16+10}" font-size="10">{BASES[bi]}</text>')
    s.append('</svg>')
    return "\n".join(s)


def fig3():
    w, h = 900, 440
    L, R, T, B = 60, 870, 70, 370
    maxd = 0.04
    s = [svg_head(w, h, "Fig. 3.  Dinucleotide deviation from independence (red band = TA under-representation)")]
    # zero line + labels
    for yv in [-0.04, -0.02, 0.0, 0.02, 0.04]:
        py = B - ((yv + maxd) / (2 * maxd)) * (B - T)
        s.append(f'<line x1="{L}" y1="{py:.1f}" x2="{R}" y2="{py:.1f}" stroke="#DDDDDD"/>')
        s.append(f'<text x="{L-6}" y="{py+4:.1f}" text-anchor="end" font-size="9">{yv:.2f}</text>')
    s.append(f'<text x="6" y="44" font-size="12">Observed - expected frequency</text>')
    s.append(f'<text x="{(L+R)/2}" y="{B+42}" text-anchor="middle" font-size="12">Dinucleotide</text>')
    s.append(f'<line x1="{L}" y1="{T}" x2="{L}" y2="{B}" stroke="black"/>')
    s.append(f'<line x1="{L}" y1="{B}" x2="{R}" y2="{B}" stroke="black"/>')
    d12 = [O12[i] - E12[i] for i in range(16)]
    d23 = [O23[i] - E23[i] for i in range(16)]
    d31 = [O31[i] - E31[i] for i in range(16)]
    gw = (R - L) / 16
    bw = gw / 4
    col12, col23, col31 = '#4682B4', '#FF8C00', '#2E8B57'
    for i in range(16):
        gx0 = L + i * gw
        if DINUC[i] == 'TA':
            s.append(f'<rect x="{gx0}" y="{T}" width="{gw}" height="{B-T}" '
                     f'fill="#FF0000" opacity="0.12"/>')
        for val, col, off in [(d12[i], col12, -1), (d23[i], col23, 0), (d31[i], col31, 1)]:
            bh = abs(val) / (2 * maxd) * (B - T)
            bx = gx0 + gw / 2 + off * bw
            by = B - bh if val >= 0 else B
            s.append(f'<rect x="{bx:.1f}" y="{by:.1f}" width="{bw-1:.1f}" height="{bh:.1f}" fill="{col}"/>')
        s.append(f'<text x="{gx0+gw/2}" y="{B+16}" text-anchor="middle" font-size="9">{DINUC[i]}</text>')
    # legend
    lx, ly = R - 150, T + 6
    for idx, (col, lab) in enumerate([(col12, 'positions 1-2'), (col23, 'positions 2-3'), (col31, 'positions 3-1')]):
        s.append(f'<rect x="{lx}" y="{ly+idx*16}" width="12" height="12" fill="{col}"/>')
        s.append(f'<text x="{lx+18}" y="{ly+idx*16+10}" font-size="9">{lab}</text>')
    s.append('</svg>')
    return "\n".join(s)


def main():
    out = {
        'Fig1_codon_frequency.svg': fig1(),
        'Fig2_base_composition.svg': fig2(),
        'Fig3_dinucleotide_deviation.svg': fig3(),
    }
    for name, content in out.items():
        with open(name, 'w', encoding='utf-8') as f:
            f.write(content)
        print(f"  wrote {name}  ({len(content)} bytes)")


if __name__ == '__main__':
    print("Generating figures (SVG, standard library only) ...")
    main()
    print("Done. Convert to PNG/PDF at 300 dpi with any SVG tool if required.")
