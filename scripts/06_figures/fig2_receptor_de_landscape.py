#!/usr/bin/env python3
"""
Figure 2: Receptor Differential Expression Landscape (Mouse)

Panel A: Volcano Plot - All Receptors (Wilcoxon scores on x-axis)
Panel B: Receptor Ranking Bar Plot (Top 20 by score)
Panel C: Pathway-Level Summary
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

# Load data
rec = pd.read_csv(PROJECT_DIR / "results" / "mouse" / "receptor_rankings_by_logfc.csv")

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

# Use Wilcoxon scores for x-axis
rec['score_capped'] = np.clip(rec['scores'], -200, 200)
rec['neg_log_padj'] = -np.log10(rec['pvals_adj'].clip(lower=1e-300))

# Highlight genes
highlights = ['Fgfr2','Bmpr2','Fzd2','Acvr1','Egfr','Tgfbr1','Notch1','Kdr','Pdgfra','Epha7']

# Color map for pathways
pw_colors = {
    'FGF':'#E74C3C','BMP':'#3498DB','Wnt':'#2ECC71','TGF-\u03b2':'#9B59B6',
    'PDGF':'#F39C12','VEGF':'#1ABC9C','Notch':'#E67E22','EGF':'#34495E',
    'HGF':'#16A085','IGF':'#8E44AD','Ephrin':'#D35400','Endothelin':'#C0392B',
    'TAM':'#7F8C8D','Sema':'#2980B9','Sema/VEGF':'#27AE60','Other':'#BDC3C7'
}

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# ---- Panel A: Volcano Plot (scores on x-axis) ----
ax = axes[0]
# Background (Other)
other = rec[rec['pathway'] == 'Other']
ax.scatter(other['score_capped'], other['neg_log_padj'],
           c='#D5D8DC', s=8, alpha=0.3, zorder=1)

# Colored by pathway
for pw, color in pw_colors.items():
    if pw == 'Other':
        continue
    sub = rec[rec['pathway'] == pw]
    if len(sub) > 0:
        ax.scatter(sub['score_capped'], sub['neg_log_padj'],
                   c=color, s=20, alpha=0.7, label=pw, zorder=2)

# Highlight specific genes
for gene in highlights:
    row = rec[rec['names'] == gene]
    if len(row) > 0:
        r = row.iloc[0]
        ax.annotate(gene, (r['score_capped'], r['neg_log_padj']),
                    fontsize=7, fontweight='bold',
                    xytext=(5, 5), textcoords='offset points',
                    arrowprops=dict(arrowstyle='-', color='black', lw=0.5),
                    zorder=3)

ax.axhline(-np.log10(0.05), color='gray', linestyle='--', linewidth=0.5)
ax.axvline(0, color='gray', linestyle='--', linewidth=0.5)
ax.set_xlabel('Wilcoxon Score (Activated vs Quiescent)', fontsize=10)
ax.set_ylabel('-log\u2081\u2080(padj)', fontsize=10)
ax.set_title('A. Receptor Volcano Plot', fontsize=12, fontweight='bold')
ax.legend(fontsize=6, loc='upper left', ncol=2, framealpha=0.8)

# ---- Panel B: Top 20 Receptors Bar Plot (by score) ----
ax = axes[1]
top20 = rec[rec['scores'] > 0].nlargest(20, 'scores').iloc[::-1]
colors = [pw_colors.get(top20.iloc[i]['pathway'], '#BDC3C7') for i in range(len(top20))]
bars = ax.barh(range(len(top20)), np.clip(top20['scores'], 0, 200), color=colors)

# Mark FGFR2 position
for i, (_, row) in enumerate(top20.iterrows()):
    ax.text(np.clip(row['scores'], 0, 200) + 1, i,
            f"{row['pathway']}", fontsize=6, va='center', color='#555')
    if row['names'] == 'Fgfr2':
        bars[i].set_edgecolor('red')
        bars[i].set_linewidth(2)

ax.set_yticks(range(len(top20)))
ax.set_yticklabels(top20['names'], fontsize=8)
ax.set_xlabel('Wilcoxon Score', fontsize=10)
ax.set_title('B. Top 20 Upregulated Receptors (by Score)', fontsize=12, fontweight='bold')

# ---- Panel C: Pathway-Level Summary ----
ax = axes[2]
sig = rec[rec['pvals_adj'] < 0.05]
pathways_of_interest = ['FGF','BMP','Wnt','TGF-\u03b2','PDGF','VEGF','Notch','EGF','Ephrin','TAM','Sema']

pw_data = []
for pw in pathways_of_interest:
    sub = sig[sig['pathway'] == pw]
    n_up = (sub['scores'] > 0).sum()
    n_down = (sub['scores'] < 0).sum()
    pw_data.append({'pathway': pw, 'up': n_up, 'down': -n_down})

pw_df = pd.DataFrame(pw_data)
x = range(len(pw_df))
ax.barh(x, pw_df['up'], color='#E74C3C', alpha=0.7, label='Upregulated')
ax.barh(x, pw_df['down'], color='#3498DB', alpha=0.7, label='Downregulated')
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
