#!/usr/bin/env python3
"""
Figure 2: Receptor Differential Expression Landscape (Mouse) — OR-based V2

Panel A: Volcano Plot - log(OR) vs -log10(Fisher p)
Panel B: Waterfall of significantly upregulated receptors ranked by log(OR)
Panel C: Pathway-Level Summary (OR > 1 up, OR < 1 down, Fisher p < 0.05)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent.parent.parent
OUTPUT_DIR = PROJECT_DIR / "results" / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Load receptor rankings (has logFC, scores, padj)
rec = pd.read_csv(PROJECT_DIR / "results" / "mouse" / "receptor_rankings_by_logfc.csv")

# Load OR data for Panel B waterfall
or_data = pd.read_csv(PROJECT_DIR / "results" / "mismatch" / "mouse_gene_odds_ratios.csv")
or_map = or_data.set_index('gene')['log_or'].to_dict()
fisher_map = or_data.set_index('gene')['fisher_p'].to_dict()
rec['log_or'] = rec['names'].map(or_map)
rec['fisher_p'] = rec['names'].map(fisher_map)

# Panel A: remove logFC artifact genes (|logFC/score| ratio > 3)
rec['logfc_score_ratio'] = rec['logfoldchanges'].abs() / (rec['scores'].abs() + 0.01)
rec['is_artifact'] = rec['logfc_score_ratio'] > 3
rec_clean = rec[~rec['is_artifact']].copy()
print(f"Filtered {rec['is_artifact'].sum()}/{len(rec)} artifact genes for volcano")

# Pathway annotation
pathway_map = {
    'Fgfr1':'FGF','Fgfr2':'FGF','Fgfr3':'FGF','Fgfrl1':'FGF',
    'Bmpr1a':'BMP','Bmpr1b':'BMP','Bmpr2':'BMP',
    'Acvr1':'BMP','Acvr1b':'BMP','Acvr2a':'BMP',
    'Tgfbr1':'TGF-\u03b2','Tgfbr2':'TGF-\u03b2','Tgfbr3':'TGF-\u03b2',
    'Fzd1':'Wnt','Fzd2':'Wnt','Fzd3':'Wnt','Fzd4':'Wnt',
    'Fzd5':'Wnt','Fzd6':'Wnt','Fzd7':'Wnt','Fzd8':'Wnt',
    'Lrp5':'Wnt','Lrp6':'Wnt','Ryk':'Wnt','Ror2':'Wnt',
    'Pdgfra':'PDGF','Pdgfrb':'PDGF','Pdgfrl':'PDGF',
    'Egfr':'EGF','Erbb2':'EGF','Erbb3':'EGF','Erbb4':'EGF',
    'Kdr':'VEGF','Flt1':'VEGF','Flt4':'VEGF',
    'Notch1':'Notch','Notch2':'Notch','Notch3':'Notch','Notch4':'Notch',
    'Met':'HGF','Igf1r':'IGF','Igf2r':'IGF',
    'Epha1':'Ephrin','Epha3':'Ephrin','Epha7':'Ephrin',
    'Ephb2':'Ephrin','Ephb3':'Ephrin','Ephb6':'Ephrin',
    'Ednra':'Endothelin','Ednrb':'Endothelin',
    'Tyro3':'TAM','Mertk':'TAM',
    'Nrp1':'Sema/VEGF','Nrp2':'Sema/VEGF',
    'Plxna4':'Sema','Plxnb2':'Sema',
}
rec['pathway'] = rec['names'].map(pathway_map).fillna('Other')
rec_clean['pathway'] = rec_clean['names'].map(pathway_map).fillna('Other')

# Volcano columns (logFC based, for Panel A)
rec_clean['logfc_capped'] = np.clip(rec_clean['logfoldchanges'], -15, 15)
rec_clean['neg_log_padj'] = -np.log10(rec_clean['pvals_adj'].clip(lower=1e-300))

# Highlights
highlights = ['Fgfr2','Bmpr2','Fzd2','Acvr1','Egfr','Tgfbr1','Notch1','Kdr','Pdgfra','Epha7']

# Color map for pathways
pw_colors = {
    'FGF':'#E74C3C','BMP':'#3498DB','Wnt':'#2ECC71','TGF-\u03b2':'#9B59B6',
    'PDGF':'#F39C12','VEGF':'#1ABC9C','Notch':'#E67E22','EGF':'#34495E',
    'HGF':'#16A085','IGF':'#8E44AD','Ephrin':'#D35400','Endothelin':'#C0392B',
    'TAM':'#7F8C8D','Sema':'#2980B9','Sema/VEGF':'#27AE60','Other':'#BDC3C7'
}

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# ---- Panel A: Volcano Plot (logFC, artifacts removed) ----
ax = axes[0]
other = rec_clean[rec_clean['pathway'] == 'Other']
ax.scatter(other['logfc_capped'], other['neg_log_padj'],
           c='#D5D8DC', s=8, alpha=0.3, zorder=1)

for pw, color in pw_colors.items():
    if pw == 'Other':
        continue
    sub = rec_clean[rec_clean['pathway'] == pw]
    if len(sub) > 0:
        ax.scatter(sub['logfc_capped'], sub['neg_log_padj'],
                   c=color, s=20, alpha=0.7, label=pw, zorder=2)

for gene in highlights:
    row = rec_clean[rec_clean['names'] == gene]
    if len(row) > 0:
        r = row.iloc[0]
        ax.annotate(gene, (r['logfc_capped'], r['neg_log_padj']),
                    fontsize=7, fontweight='bold',
                    xytext=(5, 5), textcoords='offset points',
                    arrowprops=dict(arrowstyle='-', color='black', lw=0.5),
                    zorder=3)

ax.axhline(-np.log10(0.05), color='gray', linestyle='--', linewidth=0.5)
ax.axvline(0, color='gray', linestyle='--', linewidth=0.5)
ax.set_xlabel('log₂FC (Activated vs Quiescent)', fontsize=10)
ax.set_ylabel('-log₁₀(padj)', fontsize=10)
ax.set_title('A. Receptor Volcano Plot', fontsize=12, fontweight='bold')
ax.legend(fontsize=6, loc='upper left', ncol=2, framealpha=0.8)

# ---- Panel B: Waterfall — all sig upregulated receptors by log(OR) ----
ax = axes[1]
sig_up = rec[(rec['log_or'] > 0) & (rec['fisher_p'] < 0.05)].dropna(subset=['log_or']).sort_values('log_or', ascending=False).reset_index(drop=True)
n_up = len(sig_up)

# Plot waterfall curve
ax.fill_between(range(n_up), sig_up['log_or'], alpha=0.15, color='#E74C3C')
ax.plot(range(n_up), sig_up['log_or'], color='#E74C3C', linewidth=0.8, alpha=0.6)

# Highlight key receptors
key_genes = {
    'Fgfr2': ('#E74C3C', 'FGF'),
    'Bmpr2': ('#3498DB', 'BMP'),
    'Acvr1': ('#2980B9', 'BMP'),
    'Fzd2': ('#2ECC71', 'Wnt'),
    'Egfr': ('#34495E', 'EGF'),
    'Notch1': ('#E67E22', 'Notch'),
    'Epha7': ('#D35400', 'Ephrin'),
    'Tgfbr1': ('#9B59B6', 'TGF-\u03b2'),
}

# Collect gene positions, sort by rank (left to right on curve)
gene_positions = []
for gene, (color, pw) in key_genes.items():
    idx = sig_up[sig_up['names'] == gene].index
    if len(idx) > 0:
        i = idx[0]
        val = sig_up.loc[i, 'log_or']
        gene_positions.append({'gene': gene, 'i': i, 'val': val, 'color': color, 'rank': i+1})

gene_positions.sort(key=lambda x: x['i'])

# Place labels along the curve: each label sits just above its point,
# alternating left/right to avoid overlap with neighbors.
# Since the curve descends left→right, labels naturally flow in rank order.
y_max = sig_up['log_or'].max()
min_y_gap = y_max * 0.055  # minimum vertical spacing between labels

# First pass: ideal y = point y + small offset
label_positions = []
for gp in gene_positions:
    label_positions.append({**gp, 'label_y': gp['val'] + y_max * 0.04})

# Second pass: push labels up if they overlap with the previous label
for i in range(1, len(label_positions)):
    prev_y = label_positions[i-1]['label_y']
    curr_y = label_positions[i]['label_y']
    if prev_y - curr_y < min_y_gap:
        # Current label is too close to previous; push previous up
        # Work backwards to create space
        label_positions[i-1]['label_y'] = curr_y + min_y_gap

# Third pass: final adjustment — work forward
for i in range(1, len(label_positions)):
    prev_y = label_positions[i-1]['label_y']
    curr_y = label_positions[i]['label_y']
    if prev_y - curr_y < min_y_gap:
        label_positions[i]['label_y'] = prev_y - min_y_gap

for lp in label_positions:
    ax.scatter(lp['i'], lp['val'], c=lp['color'], s=60, zorder=3,
               edgecolors='black', linewidths=0.5)
    ax.text(lp['i'] + n_up * 0.015, lp['label_y'],
            f"{lp['gene']} (#{lp['rank']})",
            fontsize=7, fontweight='bold', color=lp['color'],
            va='center', zorder=4)

ax.set_xlabel(f'Receptor Rank (1\u2013{n_up} upregulated)', fontsize=10)
ax.set_ylabel('log(OR)', fontsize=10)
ax.set_title(f'B. All {n_up} Upregulated Receptors (OR-based)', fontsize=12, fontweight='bold')
ax.set_xlim(-5, n_up + 5)

# ---- Panel C: Pathway-Level Summary (OR-based) ----
ax = axes[2]
sig = rec[rec['fisher_p'] < 0.05].dropna(subset=['log_or'])
pathways_of_interest = ['FGF','BMP','Wnt','TGF-\u03b2','PDGF','VEGF','Notch','EGF','Ephrin','TAM','Sema']

pw_data = []
for pw in pathways_of_interest:
    sub = sig[sig['pathway'] == pw]
    n_up_pw = (sub['log_or'] > 0).sum()
    n_down_pw = (sub['log_or'] < 0).sum()
    pw_data.append({'pathway': pw, 'up': n_up_pw, 'down': -n_down_pw})

pw_df = pd.DataFrame(pw_data)
x = range(len(pw_df))
ax.barh(x, pw_df['up'], color='#E74C3C', alpha=0.7, label='Upregulated (OR > 1)')
ax.barh(x, pw_df['down'], color='#3498DB', alpha=0.7, label='Downregulated (OR < 1)')
ax.set_yticks(x)
ax.set_yticklabels(pw_df['pathway'], fontsize=9)
ax.axvline(0, color='black', linewidth=0.5)
ax.set_xlabel('Number of Receptors', fontsize=10)
ax.set_title('C. Pathway-Level Receptor Changes', fontsize=12, fontweight='bold')
ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'fig2_receptor_de_landscape.png', dpi=300, bbox_inches='tight')
plt.savefig(OUTPUT_DIR / 'fig2_receptor_de_landscape.pdf', bbox_inches='tight')
print(f"Saved: fig2_receptor_de_landscape.png/pdf")
