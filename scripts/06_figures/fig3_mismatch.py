#!/usr/bin/env python3
"""
Figure 3: "Primed But Starved" Ligand-Receptor Mismatch — OR-based V2

Panel A: Concept Diagram (schematic, unchanged)
Panel B: Heatmap - Top Mismatch Pairs (mouse, log(OR) for receptor & ligand)
Panel C: Scatter Plot - Receptor log(OR) vs Ligand log(OR) (all HC pairs)
Panel D: Mismatch Score Ranking (top 20 pairs, OR composite)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent.parent.parent
OUTPUT_DIR = PROJECT_DIR / "results" / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Load OR data for mouse genes
mouse_or = pd.read_csv(PROJECT_DIR / "results" / "mismatch" / "mouse_gene_odds_ratios.csv")
mouse_or_map = mouse_or.set_index('gene')['log_or'].to_dict()

# Load 77 mismatch pairs (for pair list / pathway info)
mismatch = pd.read_csv(PROJECT_DIR / "results" / "mismatch" / "mouse_lr_mismatch_scores.csv")

# Load cross-species pairs for scatter
cross = pd.read_csv(PROJECT_DIR / "results" / "mismatch" / "cross_species_lr_mismatch_scores.csv")

# Pathway colors
pw_colors = {
    'FGF':'#E74C3C','BMP':'#3498DB','BMP/Activin':'#2980B9','Wnt':'#2ECC71',
    'TGF-\u03b2':'#9B59B6','EGF':'#34495E','Notch':'#E67E22','Sema':'#1ABC9C',
    'Sema/VEGF':'#27AE60','IGF':'#8E44AD','Other':'#95A5A6'
}

fig, axes = plt.subplots(2, 2, figsize=(16, 14))

# ---- Panel A: Concept Diagram (unchanged) ----
ax = axes[0, 0]
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)
ax.axis('off')
ax.set_title('A. "Primed But Starved" Hypothesis', fontsize=12, fontweight='bold')

# Normal
ax.text(2.5, 9.2, 'NORMAL', fontsize=11, ha='center', fontweight='bold', color='#2C3E50')
ax.text(2.5, 8.2, 'Ligand \u25cf\u25cf\u25cf\u25cf\u25cf', fontsize=9, ha='center', family='monospace', color='#27AE60')
ax.text(2.5, 7.4, '\u2193', fontsize=14, ha='center', color='gray')
ax.text(2.5, 6.6, 'Receptor \u25a2\u25a2\u25a2', fontsize=9, ha='center', family='monospace', color='#3498DB')
ax.text(2.5, 5.8, '\u2193', fontsize=14, ha='center', color='gray')
ax.text(2.5, 5.0, 'Signal \u2588\u2588\u2588\u2588', fontsize=9, ha='center', family='monospace', color='#2ECC71')
ax.text(2.5, 4.0, 'Quiescent', fontsize=10, ha='center', style='italic', color='#7F8C8D')

# Post-MI
ax.text(7.5, 9.2, 'POST-MI', fontsize=11, ha='center', fontweight='bold', color='#E74C3C')
ax.text(7.5, 8.2, 'Ligand \u25cf\u25cf', fontsize=9, ha='center', family='monospace', color='#E74C3C')
ax.text(7.5, 7.8, '(depleted)', fontsize=7, ha='center', color='#E74C3C')
ax.text(7.5, 7.0, '\u2193', fontsize=14, ha='center', color='gray')
ax.text(7.5, 6.2, 'Receptor \u25a2\u25a2\u25a2\u25a2\u25a2\u25a2\u25a2', fontsize=9, ha='center', family='monospace', color='#3498DB')
ax.text(7.5, 5.8, '(upregulated)', fontsize=7, ha='center', color='#3498DB')
ax.text(7.5, 5.0, '\u2193', fontsize=14, ha='center', color='gray')
ax.text(7.5, 4.2, 'Signal \u2588\u2588', fontsize=9, ha='center', family='monospace', color='#F39C12')
ax.text(7.5, 3.6, '(insufficient)', fontsize=7, ha='center', color='#F39C12')

# Therapeutic strategy
ax.annotate('', xy=(7.5, 2.5), xytext=(2.5, 2.5),
            arrowprops=dict(arrowstyle='->', color='#8E44AD', lw=2))
ax.text(5.0, 2.8, 'THERAPEUTIC STRATEGY', fontsize=8, ha='center', fontweight='bold', color='#8E44AD')
ax.text(5.0, 1.8, 'Deliver depleted ligands\n\u2192 engage upregulated receptors\n\u2192 enhanced activation',
        fontsize=8, ha='center', color='#8E44AD')

# Example
ax.text(5.0, 0.6, 'Example: FGF10 (\u2193) + FGFR2 (\u2191) \u2192 Deliver FGF10',
        fontsize=8, ha='center', fontweight='bold', color='#2C3E50',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#FDEBD0', edgecolor='#E67E22'))
ax.text(5.0, -0.2, '(Conceptual diagram)', fontsize=7, ha='center', style='italic', color='#7F8C8D')

# ---- Panel B: Heatmap - Top Pairs (using log(OR)) ----
ax = axes[0, 1]

# Attach log(OR) values from mouse_gene_odds_ratios to the 77 pairs
mismatch['r_log_or'] = mismatch['receptor'].map(mouse_or_map)
mismatch['l_log_or'] = mismatch['ligand'].map(mouse_or_map)

# Compute OR-based composite: (log_or_receptor - log_or_ligand) * db_weight
# db_weight = log2(n_db + 1) as a simple weight
mismatch['db_weight'] = mismatch['n_db'] / 5.0
mismatch['composite_or'] = (mismatch['r_log_or'] - mismatch['l_log_or']) * mismatch['db_weight']

# Drop pairs with missing OR
mismatch_valid = mismatch.dropna(subset=['r_log_or', 'l_log_or']).copy()
mismatch_valid = mismatch_valid.sort_values('composite_or', ascending=False).reset_index(drop=True)

# Select top 25 pairs
top_pairs = mismatch_valid.head(25).copy()
top_pairs['pair'] = top_pairs['receptor'] + '/' + top_pairs['ligand']

# Also add FGFR2/FGF10 if not in top 25
fgf10_pair = mismatch_valid[(mismatch_valid['receptor']=='Fgfr2') & (mismatch_valid['ligand']=='Fgf10')]
if len(fgf10_pair) > 0 and 'Fgfr2/Fgf10' not in top_pairs['pair'].values:
    fgf10_pair = fgf10_pair.copy()
    fgf10_pair['pair'] = 'Fgfr2/Fgf10'
    top_pairs = pd.concat([top_pairs, fgf10_pair])

heatmap_data = top_pairs[['pair', 'r_log_or', 'l_log_or']].set_index('pair')
heatmap_data.columns = ['Receptor\nlog(OR)', 'Ligand\nlog(OR)']

sns.heatmap(heatmap_data,
            cmap='RdBu_r', center=0, vmin=-2.5, vmax=2.5,
            annot=True, fmt='.2f', annot_kws={'fontsize': 7},
            ax=ax, cbar_kws={'shrink': 0.6, 'label': 'log(OR)'})
ax.set_title('B. Top Mismatch Pairs (Mouse, OR-based)', fontsize=12, fontweight='bold')
ax.tick_params(axis='y', labelsize=6)

# ---- Panel C: Scatter - Receptor log(OR) vs Ligand log(OR) ----
ax = axes[1, 0]

# Use all HC pairs from cross_species file, match gene names to mouse OR data
# Cross-species file uses human gene names (uppercase); mouse OR uses mouse names
# We need to match: cross receptor (FGFR2) -> mouse gene (Fgfr2)
# Build a case-insensitive lookup from mouse OR
mouse_or_upper = {g.upper(): v for g, v in mouse_or_map.items()}

cross['r_log_or'] = cross['receptor'].str.upper().map(mouse_or_upper)
cross['l_log_or'] = cross['ligand'].str.upper().map(mouse_or_upper)

scatter_data = cross.dropna(subset=['r_log_or', 'l_log_or']).copy()
SCORE_LIM = 2.5

# Clip for display
scatter_data = scatter_data[
    (scatter_data['r_log_or'].abs() <= SCORE_LIM * 1.5) &
    (scatter_data['l_log_or'].abs() <= SCORE_LIM * 1.5)
].copy()

ax.scatter(scatter_data['r_log_or'], scatter_data['l_log_or'],
           c='#BDC3C7', s=10, alpha=0.3, zorder=1)

# Highlight conserved pairs
conserved = scatter_data[scatter_data['conservation'].str.startswith('CONSERVED', na=False)]
ax.scatter(conserved['r_log_or'], conserved['l_log_or'],
           c='#E74C3C', s=25, alpha=0.6, zorder=2, label='Conserved')

# Highlight key pairs
key_pairs = [
    ('FGFR2',  'FGF10', (10, -12)),
    ('FGFR2',  'FGF16', (-65, -30)),
    ('BMPR2',  'BMP6',  (8, 10)),
    ('ACVR1',  'BMP6',  (-55, -12)),
    ('BMPR2',  'BMP4',  (-55, 8)),
]
for rec_name, lig_name, offset in key_pairs:
    row = scatter_data[(scatter_data['receptor']==rec_name) & (scatter_data['ligand']==lig_name)]
    if len(row) > 0:
        r = row.iloc[0]
        ax.annotate(f'{rec_name}/{lig_name}',
                    (r['r_log_or'], r['l_log_or']),
                    fontsize=6, fontweight='bold',
                    xytext=offset, textcoords='offset points',
                    arrowprops=dict(arrowstyle='->', color='black', lw=0.5),
                    zorder=3)

# Quadrant labels
ax.axhline(0, color='gray', linewidth=0.5, linestyle='--')
ax.axvline(0, color='gray', linewidth=0.5, linestyle='--')
ax.text(SCORE_LIM*0.6, SCORE_LIM*0.6, 'Both up\n(balanced)', fontsize=7, ha='center', color='#7F8C8D')
ax.text(-SCORE_LIM*0.6, SCORE_LIM*0.6, 'R\u2193 L\u2191\n(inverse)', fontsize=7, ha='center', color='#7F8C8D')
ax.text(-SCORE_LIM*0.6, -SCORE_LIM*0.6, 'Both down\n(shutdown)', fontsize=7, ha='center', color='#7F8C8D')

# Q4 highlight
ax.fill_between([0, SCORE_LIM], -SCORE_LIM, 0, alpha=0.08, color='red')
ax.text(SCORE_LIM*0.6, -SCORE_LIM*0.6, 'R\u2191 L\u2193\nPRIMED BUT\nSTARVED \u2605',
        fontsize=8, ha='center', fontweight='bold', color='#E74C3C')

ax.set_xlabel('Receptor log(OR)', fontsize=10)
ax.set_ylabel('Ligand log(OR)', fontsize=10)
ax.set_title(f'C. Receptor vs Ligand log(OR)', fontsize=12, fontweight='bold')
ax.legend(fontsize=8)
ax.set_xlim(-SCORE_LIM, SCORE_LIM)
ax.set_ylim(-SCORE_LIM, SCORE_LIM)

# ---- Panel D: Mismatch Ranking Bar Chart (OR composite) ----
ax = axes[1, 1]

# Top 20 from OR-based ranking
top20 = mismatch_valid.head(20).copy()

# Add FGFR2/FGF10 if not in top 20
fgf10 = mismatch_valid[(mismatch_valid['receptor']=='Fgfr2') & (mismatch_valid['ligand']=='Fgf10')]
if len(fgf10) > 0:
    fgf10_rank = fgf10.index[0]
    if fgf10_rank >= 20:
        top20 = pd.concat([top20, fgf10])

top20 = top20.iloc[::-1]  # reverse for horizontal bar
top20['pair'] = top20['receptor'] + '/' + top20['ligand']
colors = [pw_colors.get(top20.iloc[i]['pathway'], '#95A5A6') for i in range(len(top20))]

bars = ax.barh(range(len(top20)), top20['composite_or'], color=colors, edgecolor='white', linewidth=0.5)

# Highlight FGFR2/FGF10
for i, (_, row) in enumerate(top20.iterrows()):
    if row['receptor'] == 'Fgfr2' and row['ligand'] == 'Fgf10':
        bars[i].set_edgecolor('#FF0000')
        bars[i].set_linewidth(3)
        bars[i].set_hatch('//')
        # Add star marker at end of bar
        ax.plot(top20.iloc[i]['composite_or'], i, marker='*', color='#FF0000',
                markersize=12, zorder=5)

ax.set_yticks(range(len(top20)))
ax.set_yticklabels(top20['pair'], fontsize=7)
ax.set_xlabel('OR Composite Mismatch Score', fontsize=10)
ax.set_title('D. Top Mismatch Pairs (OR-based)', fontsize=12, fontweight='bold')

# Legend for pathways
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=pw_colors[pw], label=pw)
                   for pw in ['FGF','BMP','BMP/Activin','Wnt','EGF','Notch','Other']
                   if pw in top20['pathway'].values or pw == 'FGF']
ax.legend(handles=legend_elements, fontsize=7, loc='lower right')

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'fig3_mismatch.png', dpi=300, bbox_inches='tight')
plt.savefig(OUTPUT_DIR / 'fig3_mismatch.pdf', bbox_inches='tight')
print("Saved: fig3_mismatch.png/pdf")
