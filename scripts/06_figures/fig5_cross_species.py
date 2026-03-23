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

# ---- Panel A: Mouse vs Human Receptor Score Scatter ----
ax = axes[0, 0]

# Get unique receptors with both species data
receptors = cross.drop_duplicates(subset='receptor')[['receptor', 'm_r_score', 'h_r_score']].dropna()
LIM_5A = 15
receptors = receptors[
    (receptors['m_r_score'].abs() <= LIM_5A) &
    (receptors['h_r_score'].abs() <= LIM_5A)
].copy()
receptors['m_r_c'] = receptors['m_r_score']
receptors['h_r_c'] = receptors['h_r_score']

ax.scatter(receptors['m_r_c'], receptors['h_r_c'],
           c='#BDC3C7', s=15, alpha=0.4, zorder=1)

# Highlight: both upregulated
both_up = receptors[(receptors['m_r_score'] > 0) & (receptors['h_r_score'] > 0)]
ax.scatter(both_up['m_r_c'], both_up['h_r_c'],
           c='#E74C3C', s=25, alpha=0.6, zorder=2, label=f'Both up (n={len(both_up)})')

# Annotate key conserved upregulated receptors with manual offsets
key_labels = [
    ('FGFR2',  (5, -12)),
    ('BMPR2',  (-45, 8)),   # outside range, won't show — skip handled below
    ('ACVR1',  (-40, 8)),
    ('NOTCH1', (5, 5)),
    ('IL2RG',  (-40, -10)),
    ('IL1RL2', (5, -12)),
    ('FGFR2',  (5, -12)),
    ('FZD1',   (-30, 8)),
]
labeled = set()
for gene, offset in key_labels:
    if gene in labeled:
        continue
    row = receptors[receptors['receptor'] == gene]
    if len(row) > 0:
        r = row.iloc[0]
        if r['m_r_score'] > 0 and r['h_r_score'] > 0:
            ax.scatter(r['m_r_c'], r['h_r_c'], c='darkred', s=50, zorder=3,
                       edgecolors='black', linewidths=0.5)
            ax.annotate(gene, (r['m_r_c'], r['h_r_c']),
                        fontsize=7, fontweight='bold', color='darkred',
                        xytext=offset, textcoords='offset points',
                        arrowprops=dict(arrowstyle='->', color='darkred', lw=0.8),
                        zorder=4)
            labeled.add(gene)

# Highlight Q1 (both upregulated) region
ax.fill_between([0, LIM_5A], 0, LIM_5A, alpha=0.05, color='red')
ax.text(LIM_5A*0.7, LIM_5A*0.8, 'Conserved\nupregulated', fontsize=7,
        ha='center', color='#C0392B', fontstyle='italic')

ax.axhline(0, color='gray', linewidth=0.5, linestyle='--')
ax.axvline(0, color='gray', linewidth=0.5, linestyle='--')
ax.set_xlabel('Mouse Receptor Score', fontsize=10)
ax.set_ylabel('Human Receptor Score', fontsize=10)
ax.set_title('A. Receptor Score: Mouse vs Human (|score| ≤ 15)', fontsize=12, fontweight='bold')
ax.legend(fontsize=8)
ax.set_xlim(-LIM_5A, LIM_5A)
ax.set_ylim(-LIM_5A, LIM_5A)

# ---- Panel B: Conservation Heatmap (using scores) ----
ax = axes[0, 1]

# Top 20 conserved pairs (by avg mismatch score)
conserved = cross[cross['conservation'].str.startswith('CONSERVED', na=False)].copy()
conserved['avg_mm'] = conserved['avg_mismatch_score'].astype(float)
conserved = conserved.sort_values('avg_mm', ascending=False)

top20 = conserved.head(20).copy()
top20['pair'] = top20['receptor'] + '/' + top20['ligand']

# Build heatmap data using scores
hm_data = pd.DataFrame({
    'Mouse R': np.clip(top20['m_r_score'].astype(float), -100, 100).values,
    'Mouse L': np.clip(top20['m_l_score'].astype(float), -100, 100).values,
    'Human R': np.clip(top20['h_r_score'].astype(float), -20, 20).values,
    'Human L': np.clip(top20['h_l_score'].astype(float), -20, 20).values,
}, index=top20['pair'].values)

sns.heatmap(hm_data, cmap='RdBu_r', center=0, vmin=-100, vmax=100,
            annot=True, fmt='.1f', annot_kws={'fontsize': 7},
            ax=ax, cbar_kws={'shrink': 0.6, 'label': 'Wilcoxon Score'})
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
