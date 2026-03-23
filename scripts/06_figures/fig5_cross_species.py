#!/usr/bin/env python3
"""
Figure 5: Cross-Species Conservation

Panel A: Mouse vs Human Score Scatter (receptors)
Panel B: Conservation Heatmap (top 20 pairs, Wilcoxon scores)
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

# Load score-based data
cross = pd.read_csv(PROJECT_DIR / "results" / "mismatch" / "cross_species_lr_mismatch_scores.csv")
targets = pd.read_csv(PROJECT_DIR / "results" / "mismatch" / "therapeutic_targets_scores.csv")

fig, axes = plt.subplots(2, 2, figsize=(16, 14))

# ---- Panel A: Mouse vs Human Receptor Effect Size (rank-biserial r) ----
ax = axes[0, 0]

# Convert z-scores to rank-biserial correlation: r = z / sqrt(N)
# This removes sample-size dependence, making cross-species comparison fair.
# Mouse: 49,366 quiescent + 63,310 activated = 112,676
# Human: 18,091 quiescent + 1,916 activated = 20,007
MOUSE_N = 112676
HUMAN_N = 20007

receptors = cross.drop_duplicates(subset='receptor')[['receptor', 'm_r_score', 'h_r_score']].dropna()
receptors['m_r_r'] = receptors['m_r_score'] / np.sqrt(MOUSE_N)
receptors['h_r_r'] = receptors['h_r_score'] / np.sqrt(HUMAN_N)

ax.scatter(receptors['m_r_r'], receptors['h_r_r'],
           c='#BDC3C7', s=15, alpha=0.4, zorder=1)

# Highlight: both upregulated
both_up = receptors[(receptors['m_r_r'] > 0) & (receptors['h_r_r'] > 0)]
ax.scatter(both_up['m_r_r'], both_up['h_r_r'],
           c='#E74C3C', s=25, alpha=0.6, zorder=2, label=f'Both up (n={len(both_up)})')

# Annotate key conserved upregulated receptors
key_labels = [
    ('FGFR2',  (8, -12)),
    ('BMPR2',  (8, 5)),
    ('ACVR1',  (8, -10)),
    ('BMPR1A', (-55, -12)),
    ('NOTCH1', (8, 5)),
    ('LRP6',   (-40, 8)),
    ('IL2RG',  (-40, -8)),
    ('IL1RL2', (8, -12)),
    ('FZD1',   (-30, 8)),
]
labeled = set()
for gene, offset in key_labels:
    if gene in labeled:
        continue
    row = receptors[receptors['receptor'] == gene]
    if len(row) > 0:
        r = row.iloc[0]
        if r['m_r_r'] > 0 and r['h_r_r'] > 0:
            ax.scatter(r['m_r_r'], r['h_r_r'], c='darkred', s=50, zorder=3,
                       edgecolors='black', linewidths=0.5)
            ax.annotate(gene, (r['m_r_r'], r['h_r_r']),
                        fontsize=7, fontweight='bold', color='darkred',
                        xytext=offset, textcoords='offset points',
                        arrowprops=dict(arrowstyle='->', color='darkred', lw=0.8),
                        zorder=4)
            labeled.add(gene)

# Highlight Q1
r_max = max(receptors['m_r_r'].max(), receptors['h_r_r'].max()) * 1.1
ax.fill_between([0, r_max], 0, r_max, alpha=0.05, color='red')

ax.axhline(0, color='gray', linewidth=0.5, linestyle='--')
ax.axvline(0, color='gray', linewidth=0.5, linestyle='--')

# Add diagonal reference line (equal effect size)
lim = max(abs(receptors['m_r_r']).max(), abs(receptors['h_r_r']).max()) * 1.05
ax.plot([-lim, lim], [-lim, lim], 'k--', linewidth=0.5, alpha=0.3)

ax.set_xlim(-0.2, 0.2)
ax.set_ylim(-0.1, 0.1)
ax.set_xlabel('Mouse Effect Size (rank-biserial r)', fontsize=10)
ax.set_ylabel('Human Effect Size (rank-biserial r)', fontsize=10)
ax.set_title('A. Receptor Effect Size: Mouse vs Human', fontsize=12, fontweight='bold')
ax.legend(fontsize=8)

# ---- Panel B: Conservation Heatmap (using scores) ----
ax = axes[0, 1]

# Top 20 conserved pairs (by avg mismatch score)
conserved = cross[cross['conservation'].str.startswith('CONSERVED', na=False)].copy()
conserved['avg_mm'] = conserved['avg_mismatch_score'].astype(float)
conserved = conserved.sort_values('avg_mm', ascending=False)

top20 = conserved.head(20).copy()
top20['pair'] = top20['receptor'] + '/' + top20['ligand']

# Build heatmap data using rank-biserial r = z / sqrt(N)
hm_data = pd.DataFrame({
    'Mouse R': (top20['m_r_score'].astype(float) / np.sqrt(MOUSE_N)).values,
    'Mouse L': (top20['m_l_score'].astype(float) / np.sqrt(MOUSE_N)).values,
    'Human R': (top20['h_r_score'].astype(float) / np.sqrt(HUMAN_N)).values,
    'Human L': (top20['h_l_score'].astype(float) / np.sqrt(HUMAN_N)).values,
}, index=top20['pair'].values)

sns.heatmap(hm_data, cmap='RdBu_r', center=0, vmin=-0.3, vmax=0.3,
            annot=True, fmt='.3f', annot_kws={'fontsize': 6},
            ax=ax, cbar_kws={'shrink': 0.6, 'label': 'Effect Size (rank-biserial r)'})
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
ax.set_title('D. Top 10 Therapeutic Targets (Score-Based)', fontsize=12, fontweight='bold')

top10 = targets.head(10)
table_data = []
for i, (_, r) in enumerate(top10.iterrows()):
    table_data.append([
        i + 1,
        f"{r['receptor']}/{r['ligand']}",
        f"{r['priority_score']:.3f}",
        r['conservation'].replace('CONSERVED_', '').replace('mouse_only', 'mouse'),
        f"{r['druggability']:.2f}",
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
