#!/usr/bin/env python3
"""
Figure 5: Cross-Species Conservation

Panel A: Mouse vs Human logFC Scatter (receptors)
Panel B: Conservation Heatmap (top 20 pairs)
Panel C: Venn Diagram (mouse-only / human-only / conserved)
Panel D: Final Prioritized Targets Table
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
from matplotlib_venn import venn2
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent.parent.parent
OUTPUT_DIR = PROJECT_DIR / "results" / "figures"

# Load data
cross = pd.read_csv(PROJECT_DIR / "results" / "mismatch" / "cross_species_lr_mismatch.csv")
targets = pd.read_csv(PROJECT_DIR / "results" / "mismatch" / "therapeutic_targets_corrected.csv")

fig, axes = plt.subplots(2, 2, figsize=(16, 14))

# ---- Panel A: Mouse vs Human Receptor logFC Scatter ----
ax = axes[0, 0]

# Get unique receptors with both species data
receptors = cross.drop_duplicates(subset='receptor')[['receptor', 'm_r_lfc', 'h_r_lfc']].dropna()
receptors['m_r_c'] = np.clip(receptors['m_r_lfc'], -15, 15)
receptors['h_r_c'] = np.clip(receptors['h_r_lfc'], -5, 5)

ax.scatter(receptors['m_r_c'], receptors['h_r_c'],
           c='#BDC3C7', s=15, alpha=0.4, zorder=1)

# Highlight: both upregulated
both_up = receptors[(receptors['m_r_lfc'] > 0) & (receptors['h_r_lfc'] > 0)]
ax.scatter(both_up['m_r_c'], both_up['h_r_c'],
           c='#E74C3C', s=25, alpha=0.6, zorder=2, label=f'Both up (n={len(both_up)})')

# Annotate key receptors
for gene in ['FGFR2','BMPR2','ACVR1','EPHA7','TYRO3','NOTCH1','EDNRA','LRP6','INSR']:
    row = receptors[receptors['receptor'] == gene]
    if len(row) > 0:
        r = row.iloc[0]
        if r['m_r_lfc'] > 0 and r['h_r_lfc'] > 0:
            ax.annotate(gene, (r['m_r_c'], r['h_r_c']),
                        fontsize=7, fontweight='bold',
                        xytext=(5, 5), textcoords='offset points',
                        arrowprops=dict(arrowstyle='-', color='black', lw=0.5))

ax.axhline(0, color='gray', linewidth=0.5, linestyle='--')
ax.axvline(0, color='gray', linewidth=0.5, linestyle='--')
ax.set_xlabel('Mouse Receptor logFC', fontsize=10)
ax.set_ylabel('Human Receptor logFC', fontsize=10)
ax.set_title('A. Receptor logFC: Mouse vs Human', fontsize=12, fontweight='bold')
ax.legend(fontsize=8)

# ---- Panel B: Conservation Heatmap ----
ax = axes[0, 1]

# Top 20 conserved canonical pairs (by avg mismatch, family-matched)
conserved = cross[cross['conservation'].str.startswith('CONSERVED', na=False)].copy()
conserved['avg_mm'] = (
    np.clip(conserved['m_r_lfc'].fillna(0), 0, 10) - np.clip(conserved['m_l_lfc'].fillna(0), -10, 0) +
    np.clip(conserved['h_r_lfc'].fillna(0), 0, 10) - np.clip(conserved['h_l_lfc'].fillna(0), -10, 0)
) / 2
conserved = conserved.sort_values('avg_mm', ascending=False)

top20 = conserved.head(20).copy()
top20['pair'] = top20['receptor'] + '/' + top20['ligand']

# Build heatmap data
hm_data = pd.DataFrame({
    'Mouse R': np.clip(top20['m_r_lfc'], -10, 10).values,
    'Mouse L': np.clip(top20['m_l_lfc'], -10, 10).values,
    'Human R': np.clip(top20['h_r_lfc'], -5, 5).values,
    'Human L': np.clip(top20['h_l_lfc'], -5, 5).values,
}, index=top20['pair'].values)

sns.heatmap(hm_data, cmap='RdBu_r', center=0, vmin=-10, vmax=10,
            annot=True, fmt='.1f', annot_kws={'fontsize': 7},
            ax=ax, cbar_kws={'shrink': 0.6, 'label': 'logFC'})
ax.set_title('B. Conservation Heatmap (Top 20)', fontsize=12, fontweight='bold')
ax.tick_params(axis='y', labelsize=7)

# ---- Panel C: Venn Diagram ----
ax = axes[1, 0]

# Count pairs by conservation category
n_mouse_only = (cross['conservation'] == 'mouse_only').sum()
n_human_only = (cross['conservation'] == 'human_only').sum()
n_conserved = cross['conservation'].str.startswith('CONSERVED', na=False).sum()

try:
    v = venn2(subsets=(n_mouse_only, n_human_only, n_conserved),
              set_labels=('Mouse\nmismatch', 'Human\nmismatch'),
              set_colors=('#3498DB', '#E74C3C'), alpha=0.6, ax=ax)
    if v.get_label_by_id('10'):
        v.get_label_by_id('10').set_text(f'{n_mouse_only}')
    if v.get_label_by_id('01'):
        v.get_label_by_id('01').set_text(f'{n_human_only}')
    if v.get_label_by_id('11'):
        v.get_label_by_id('11').set_text(f'{n_conserved}\nConserved')
except Exception:
    # Fallback if matplotlib_venn not available
    ax.text(0.5, 0.5, f'Mouse only: {n_mouse_only}\nHuman only: {n_human_only}\nConserved: {n_conserved}',
            transform=ax.transAxes, fontsize=12, ha='center', va='center')

ax.set_title('C. Cross-Species Mismatch Overlap', fontsize=12, fontweight='bold')

# ---- Panel D: Final Targets Table ----
ax = axes[1, 1]
ax.axis('off')
ax.set_title('D. Top 10 Therapeutic Targets (Automated)', fontsize=12, fontweight='bold')

top10 = targets.head(10)
table_data = []
for i, (_, r) in enumerate(top10.iterrows()):
    table_data.append([
        i + 1,
        f"{r['receptor']}/{r['ligand']}",
        f"{r['priority_score']:.3f}",
        r['conservation'].replace('CONSERVED_', '').replace('mouse_only', 'mouse'),
        f"{r['druggability_corrected']:.2f}",
        r.get('pathway', ''),
    ])

table = ax.table(
    cellText=table_data,
    colLabels=['Rank', 'Receptor/Ligand', 'Score', 'Conservation', 'Druggability', 'Pathway'],
    loc='center',
    cellLoc='center',
)
table.auto_set_font_size(False)
table.set_fontsize(8)
table.auto_set_column_width([0, 1, 2, 3, 4, 5])
table.scale(1, 1.5)

# Color header
for j in range(6):
    table[0, j].set_facecolor('#2C3E50')
    table[0, j].set_text_props(color='white', fontweight='bold')

# Highlight FGFR2 rows
for i, (_, r) in enumerate(top10.iterrows()):
    if 'Fgfr2' in str(r['receptor']):
        for j in range(6):
            table[i + 1, j].set_facecolor('#FDEBD0')

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'fig5_cross_species.png', dpi=300, bbox_inches='tight')
plt.savefig(OUTPUT_DIR / 'fig5_cross_species.pdf', bbox_inches='tight')
print("Saved: fig5_cross_species.png/pdf")
