#!/usr/bin/env python3
"""
Figure 3: "Primed But Starved" Ligand-Receptor Mismatch

Panel A: Concept Diagram (schematic)
Panel B: Heatmap - Top Mismatch Pairs (mouse, Wilcoxon scores)
Panel C: Scatter Plot - Receptor vs Ligand Score (quadrant plot)
Panel D: Mismatch Score Ranking (top 20 pairs, score-based composite)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
from matplotlib.patches import FancyArrowPatch
import seaborn as sns
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent.parent.parent
OUTPUT_DIR = PROJECT_DIR / "results" / "figures"

# Load score-based data
mismatch = pd.read_csv(PROJECT_DIR / "results" / "mismatch" / "mouse_lr_mismatch_scores.csv")
cross = pd.read_csv(PROJECT_DIR / "results" / "mismatch" / "cross_species_lr_mismatch_scores.csv")

# Load full mouse DEG for the scatter plot
deg = pd.read_csv(PROJECT_DIR / "results" / "mouse" / "quaife_ryan_de_activated_vs_quiescent.csv")
deg = deg.dropna(subset=['scores'])

# Pathway colors
pw_colors = {
    'FGF':'#E74C3C','BMP':'#3498DB','BMP/Activin':'#2980B9','Wnt':'#2ECC71',
    'TGF-\u03b2':'#9B59B6','EGF':'#34495E','Notch':'#E67E22','Sema':'#1ABC9C',
    'Sema/VEGF':'#27AE60','IGF':'#8E44AD','Other':'#95A5A6'
}

fig, axes = plt.subplots(2, 2, figsize=(16, 14))

# ---- Panel A: Concept Diagram ----
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

# ---- Panel B: Heatmap - Top Pairs (using Wilcoxon scores) ----
ax = axes[0, 1]

# Select top 25 pairs from score-based analysis
top_pairs = mismatch.head(25).copy()
top_pairs['pair'] = top_pairs['receptor'] + '/' + top_pairs['ligand']
top_pairs['r_score_c'] = np.clip(top_pairs['r_score'], -100, 100)
top_pairs['l_score_c'] = np.clip(top_pairs['l_score'], -100, 100)

# Also add FGFR2/FGF10 if not in top 25
fgf10_pair = mismatch[(mismatch['receptor']=='Fgfr2') & (mismatch['ligand']=='Fgf10')]
if len(fgf10_pair) > 0 and 'Fgfr2/Fgf10' not in top_pairs['pair'].values:
    fgf10_pair = fgf10_pair.copy()
    fgf10_pair['pair'] = 'Fgfr2/Fgf10'
    fgf10_pair['r_score_c'] = np.clip(fgf10_pair['r_score'], -100, 100)
    fgf10_pair['l_score_c'] = np.clip(fgf10_pair['l_score'], -100, 100)
    top_pairs = pd.concat([top_pairs, fgf10_pair])

heatmap_data = top_pairs[['pair', 'r_score_c', 'l_score_c', 'composite']].set_index('pair')
heatmap_data.columns = ['Receptor\nScore', 'Ligand\nScore', 'Composite\nScore']

sns.heatmap(heatmap_data[['Receptor\nScore', 'Ligand\nScore']],
            cmap='RdBu_r', center=0, vmin=-100, vmax=100,
            annot=True, fmt='.1f', annot_kws={'fontsize': 7},
            ax=ax, cbar_kws={'shrink': 0.6, 'label': 'Wilcoxon Score'})
ax.set_title('B. Top Mismatch Pairs (Mouse)', fontsize=12, fontweight='bold')
ax.tick_params(axis='y', labelsize=7)

# ---- Panel C: Scatter - Receptor vs Ligand Score ----
ax = axes[1, 0]

# Get all L-R pairs from cross-species data (mouse score values)
scatter_data = cross[cross['m_r_score'].notna() & cross['m_l_score'].notna()].copy()
scatter_data['m_r_c'] = np.clip(scatter_data['m_r_score'], -200, 200)
scatter_data['m_l_c'] = np.clip(scatter_data['m_l_score'], -200, 200)

ax.scatter(scatter_data['m_r_c'], scatter_data['m_l_c'],
           c='#BDC3C7', s=10, alpha=0.3, zorder=1)

# Highlight conserved pairs
conserved = scatter_data[scatter_data['conservation'].str.startswith('CONSERVED', na=False)]
ax.scatter(conserved['m_r_c'], conserved['m_l_c'],
           c='#E74C3C', s=25, alpha=0.6, zorder=2, label='Conserved')

# Highlight key pairs
key_pairs = [('FGFR2','FGF10'),('BMPR2','BMP6'),('ACVR1','BMP6'),('TNFRSF12A','TNFSF12'),('ITGB1','LAMC2')]
for rec_name, lig_name in key_pairs:
    row = scatter_data[(scatter_data['receptor']==rec_name) & (scatter_data['ligand']==lig_name)]
    if len(row) > 0:
        r = row.iloc[0]
        ax.annotate(f'{rec_name}/{lig_name}',
                    (r['m_r_c'], r['m_l_c']),
                    fontsize=6, fontweight='bold',
                    xytext=(8, 8), textcoords='offset points',
                    arrowprops=dict(arrowstyle='->', color='black', lw=0.5),
                    zorder=3)

# Quadrant labels
ax.axhline(0, color='gray', linewidth=0.5, linestyle='--')
ax.axvline(0, color='gray', linewidth=0.5, linestyle='--')
ax.text(150, 150, 'Both up\n(balanced)', fontsize=7, ha='center', color='#7F8C8D')
ax.text(-150, 150, 'R\u2193 L\u2191\n(inverse)', fontsize=7, ha='center', color='#7F8C8D')
ax.text(-150, -150, 'Both down\n(shutdown)', fontsize=7, ha='center', color='#7F8C8D')

# Q4 highlight
ax.fill_between([0, 200], -200, 0, alpha=0.08, color='red')
ax.text(150, -150, 'R\u2191 L\u2193\nPRIMED BUT\nSTARVED \u2605',
        fontsize=8, ha='center', fontweight='bold', color='#E74C3C')

ax.set_xlabel('Receptor Score (Activated vs Quiescent)', fontsize=10)
ax.set_ylabel('Ligand Score (Activated vs Quiescent)', fontsize=10)
ax.set_title('C. Receptor vs Ligand Score (Mouse)', fontsize=12, fontweight='bold')
ax.legend(fontsize=8)
ax.set_xlim(-200, 200)
ax.set_ylim(-200, 200)

# ---- Panel D: Mismatch Ranking Bar Chart ----
ax = axes[1, 1]

# Top 20 from score-based analysis
top20 = mismatch.head(20).copy()

# Add FGFR2/FGF10 if not in top 20
fgf10 = mismatch[(mismatch['receptor']=='Fgfr2') & (mismatch['ligand']=='Fgf10')]
if len(fgf10) > 0:
    fgf10_rank = fgf10.index[0]
    if fgf10_rank > 20:
        top20 = pd.concat([top20, fgf10])

top20 = top20.iloc[::-1]  # reverse for horizontal bar
top20['pair'] = top20['receptor'] + '/' + top20['ligand']
colors = [pw_colors.get(top20.iloc[i]['pathway'], '#95A5A6') for i in range(len(top20))]

bars = ax.barh(range(len(top20)), top20['composite'], color=colors, edgecolor='white', linewidth=0.5)

# Highlight FGFR2/FGF10
for i, (_, row) in enumerate(top20.iterrows()):
    if row['receptor'] == 'Fgfr2' and row['ligand'] == 'Fgf10':
        bars[i].set_edgecolor('#E74C3C')
        bars[i].set_linewidth(2)
        bars[i].set_hatch('//')

ax.set_yticks(range(len(top20)))
ax.set_yticklabels(top20['pair'], fontsize=7)
ax.set_xlabel('Composite Mismatch Score (Score-Based)', fontsize=10)
ax.set_title('D. Top Mismatch Pairs (Score-Based)', fontsize=12, fontweight='bold')

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
