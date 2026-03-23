#!/usr/bin/env python3
"""
Figure 5: Cross-Species Conservation — OR-based V2

Panel A: Mouse vs Human receptor log(OR) scatter
Panel B: Conservation Heatmap (top 20 pairs, log(OR) for R and L in both species)
Panel C: Venn Diagram (unchanged)
Panel D: Final Prioritized Targets Table (OR composite scoring)
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

# Load OR data
mouse_or = pd.read_csv(PROJECT_DIR / "results" / "mismatch" / "mouse_gene_odds_ratios.csv")
human_or = pd.read_csv(PROJECT_DIR / "results" / "mismatch" / "human_gene_odds_ratios.csv")

# Build lookup dicts (mouse names like Fgfr2, human names like FGFR2)
mouse_or_map = mouse_or.set_index('gene')['log_or'].to_dict()
human_or_map = human_or.set_index('gene')['log_or'].to_dict()

# Mouse gene -> uppercase for cross-species matching
mouse_or_upper = {g.upper(): v for g, v in mouse_or_map.items()}

# Load cross-species data (for pair info, conservation labels)
cross = pd.read_csv(PROJECT_DIR / "results" / "mismatch" / "cross_species_lr_mismatch_scores.csv")

# Load mismatch pairs (for pathway info)
mismatch = pd.read_csv(PROJECT_DIR / "results" / "mismatch" / "mouse_lr_mismatch_scores.csv")

# Load gene scores for druggability/literature
gene_scores = pd.read_csv(PROJECT_DIR / "results" / "mismatch" / "gene_scores_corrected.csv")
# First column might be unnamed index — use 'gene' column as index
if 'gene' in gene_scores.columns:
    gene_scores = gene_scores.set_index('gene')
else:
    gene_scores = gene_scores.set_index(gene_scores.columns[0])
gene_scores_map = gene_scores.to_dict('index')

try:
    from matplotlib_venn import venn2
    HAS_VENN = True
except ImportError:
    HAS_VENN = False

fig, axes = plt.subplots(2, 2, figsize=(16, 14))

# ---- Panel A: Mouse vs Human Receptor log(OR) ----
ax = axes[0, 0]

# Get unique receptors from cross-species data (human gene names)
receptors = cross.drop_duplicates(subset='receptor')[['receptor']].copy()
receptors['m_r_log_or'] = receptors['receptor'].str.upper().map(mouse_or_upper)
receptors['h_r_log_or'] = receptors['receptor'].map(human_or_map)

receptors = receptors.dropna(subset=['m_r_log_or', 'h_r_log_or']).copy()
print(f"Receptors with OR data in both species: {len(receptors)}")

ax.scatter(receptors['m_r_log_or'], receptors['h_r_log_or'],
           c='#BDC3C7', s=15, alpha=0.4, zorder=1)

# Highlight: both upregulated
both_up = receptors[(receptors['m_r_log_or'] > 0) & (receptors['h_r_log_or'] > 0)]
ax.scatter(both_up['m_r_log_or'], both_up['h_r_log_or'],
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
        if r['m_r_log_or'] > 0 and r['h_r_log_or'] > 0:
            ax.scatter(r['m_r_log_or'], r['h_r_log_or'], c='darkred', s=50, zorder=3,
                       edgecolors='black', linewidths=0.5)
            ax.annotate(gene, (r['m_r_log_or'], r['h_r_log_or']),
                        fontsize=7, fontweight='bold', color='darkred',
                        xytext=offset, textcoords='offset points',
                        arrowprops=dict(arrowstyle='->', color='darkred', lw=0.8),
                        zorder=4)
            labeled.add(gene)

# Highlight Q1
axis_lim = 2.0
ax.fill_between([0, axis_lim], 0, axis_lim, alpha=0.05, color='red')

ax.axhline(0, color='gray', linewidth=0.5, linestyle='--')
ax.axvline(0, color='gray', linewidth=0.5, linestyle='--')

# Diagonal reference
ax.plot([-axis_lim, axis_lim], [-axis_lim, axis_lim], 'k--', linewidth=0.5, alpha=0.3)

ax.set_xlim(-axis_lim, axis_lim)
ax.set_ylim(-axis_lim, axis_lim)
ax.set_xlabel('Mouse Receptor log(OR)', fontsize=10)
ax.set_ylabel('Human Receptor log(OR)', fontsize=10)
ax.set_title('A. Receptor log(OR): Mouse vs Human', fontsize=12, fontweight='bold')
ax.legend(fontsize=8)

# ---- Panel B: Conservation Heatmap (log(OR), top 20 conserved pairs) ----
ax = axes[0, 1]

conserved = cross[cross['conservation'].str.startswith('CONSERVED', na=False)].copy()

# Attach log(OR) for all 4 columns
conserved['m_r_log_or'] = conserved['receptor'].str.upper().map(mouse_or_upper)
conserved['m_l_log_or'] = conserved['ligand'].str.upper().map(mouse_or_upper)
conserved['h_r_log_or'] = conserved['receptor'].map(human_or_map)
conserved['h_l_log_or'] = conserved['ligand'].map(human_or_map)

# Filter out pairs with undetected genes (log_or == -10)
conserved = conserved[
    (conserved['m_r_log_or'] > -5) & (conserved['m_l_log_or'] > -5) &
    (conserved['h_r_log_or'] > -5) & (conserved['h_l_log_or'] > -5)
].copy()

# Compute avg mismatch for sorting
conserved['avg_mm_or'] = (
    (conserved['m_r_log_or'] - conserved['m_l_log_or']) +
    (conserved['h_r_log_or'] - conserved['h_l_log_or'])
) / 2

conserved = conserved.sort_values('avg_mm_or', ascending=False)
top20 = conserved.head(20).copy()
top20['pair'] = top20['receptor'] + '/' + top20['ligand']

hm_data = pd.DataFrame({
    'Mouse R': np.clip(top20['m_r_log_or'].values, -3, 3),
    'Mouse L': np.clip(top20['m_l_log_or'].values, -3, 3),
    'Human R': np.clip(top20['h_r_log_or'].values, -3, 3),
    'Human L': np.clip(top20['h_l_log_or'].values, -3, 3),
}, index=top20['pair'].values)

sns.heatmap(hm_data, cmap='RdBu_r', center=0, vmin=-2.5, vmax=2.5,
            annot=True, fmt='.2f', annot_kws={'fontsize': 6},
            ax=ax, cbar_kws={'shrink': 0.6, 'label': 'log(OR)'})
ax.set_title('B. Conservation Heatmap (Top 20, OR-based)', fontsize=12, fontweight='bold')
ax.tick_params(axis='y', labelsize=7)

# ---- Panel C: Venn Diagram ----
ax = axes[1, 0]

n_mouse_only = (cross['conservation'] == 'mouse_only').sum()
n_human_only = (cross['conservation'] == 'human_only').sum()
n_conserved = cross['conservation'].str.startswith('CONSERVED', na=False).sum()

if HAS_VENN:
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
        ax.text(0.5, 0.5, f'Mouse only: {n_mouse_only}\nHuman only: {n_human_only}\nConserved: {n_conserved}',
                transform=ax.transAxes, fontsize=12, ha='center', va='center')
else:
    ax.text(0.5, 0.5, f'Mouse only: {n_mouse_only}\nHuman only: {n_human_only}\nConserved: {n_conserved}',
            transform=ax.transAxes, fontsize=12, ha='center', va='center')

ax.set_title('C. Cross-Species Mismatch Overlap', fontsize=12, fontweight='bold')

# ---- Panel D: Final Targets Table (OR-based priority score) ----
ax = axes[1, 1]
ax.axis('off')
ax.set_title('D. Top 10 Therapeutic Targets (OR-based)', fontsize=12, fontweight='bold')

# Recompute priority score using OR composite
# Start from the 77 mismatch pairs, attach OR values
mismatch['r_log_or'] = mismatch['receptor'].map(mouse_or_map)
mismatch['l_log_or'] = mismatch['ligand'].map(mouse_or_map)
mismatch['db_weight_or'] = mismatch['n_db'] / 5.0
mismatch['composite_or'] = (mismatch['r_log_or'] - mismatch['l_log_or']) * mismatch['db_weight_or']
mismatch_valid = mismatch.dropna(subset=['r_log_or', 'l_log_or']).copy()
mismatch_valid = mismatch_valid.sort_values('composite_or', ascending=False).reset_index(drop=True)

scored = mismatch_valid.copy()
scored['h_receptor'] = scored['receptor'].str.upper()
scored['h_ligand'] = scored['ligand'].str.upper()

# OR composite already computed: composite_or = (r_log_or - l_log_or) * db_weight
# Normalize to [0, 1]
or_max = scored['composite_or'].max()
or_min = scored['composite_or'].min()
scored['or_norm'] = (scored['composite_or'] - or_min) / (or_max - or_min) if or_max != or_min else 0

# Conservation from cross-species data
cons_map = {}
for _, row in cross.iterrows():
    cons = str(row['conservation'])
    if cons.startswith('CONSERVED_strict'):
        cons_map[(row['receptor'], row['ligand'])] = 1.0
    elif cons.startswith('CONSERVED_relaxed'):
        cons_map[(row['receptor'], row['ligand'])] = 0.7
    elif cons == 'mouse_only':
        cons_map[(row['receptor'], row['ligand'])] = 0.2
    else:
        cons_map[(row['receptor'], row['ligand'])] = 0.1

scored['cons_score'] = scored.apply(
    lambda r: cons_map.get((r['h_receptor'], r['h_ligand']), 0.1), axis=1)

# Also attach conservation label for the table
cons_label_map = {}
for _, row in cross.iterrows():
    cons_label_map[(row['receptor'], row['ligand'])] = str(row['conservation'])
scored['conservation'] = scored.apply(
    lambda r: cons_label_map.get((r['h_receptor'], r['h_ligand']), 'not_checked'), axis=1)

cross_scored = scored  # alias for downstream code

# Druggability and literature from gene_scores_corrected.csv
def get_gene_score(gene, field, default=0):
    info = gene_scores_map.get(gene, None)
    if info is None:
        return default
    val = info.get(field, default)
    return val if pd.notna(val) else default

cross_scored['druggability'] = cross_scored.apply(
    lambda r: 0.3 * get_gene_score(r['h_receptor'], 'druggability_auto', 0.1) +
              0.7 * (0.6 * get_gene_score(r['h_ligand'], 'deliver_score', 0.3) +
                     0.4 * get_gene_score(r['h_ligand'], 'druggability_auto', 0.1)),
    axis=1)

cross_scored['literature'] = cross_scored.apply(
    lambda r: (get_gene_score(r['h_receptor'], 'literature_auto', 0) +
               get_gene_score(r['h_ligand'], 'literature_auto', 0)) / 2,
    axis=1)

# Normalize druggability to [0, 1]
drug_max = cross_scored['druggability'].max()
if drug_max > 0:
    cross_scored['drug_norm'] = cross_scored['druggability'] / drug_max
else:
    cross_scored['drug_norm'] = 0

# Priority score: OR composite (30%) + conservation (30%) + druggability (20%) + literature (20%)
cross_scored['priority_score'] = (
    0.3 * cross_scored['or_norm'] +
    0.3 * cross_scored['cons_score'] +
    0.2 * cross_scored['drug_norm'] +
    0.2 * cross_scored['literature']
)

# All pairs already have positive OR composite (filtered in mismatch step)
candidates = cross_scored.sort_values('priority_score', ascending=False)
top10 = candidates.head(10)

# Pathway from mismatch data
mismatch_pw = dict(zip(
    zip(mismatch['receptor'], mismatch['ligand']),
    mismatch['pathway']
))

table_data = []
for i, (_, r) in enumerate(top10.iterrows()):
    cons_label = str(r['conservation']).replace('CONSERVED_', '').replace('mouse_only', 'mouse')
    # Try to get pathway: match receptor name (uppercase to mouse format)
    pw = ''
    for (mr, ml), p in mismatch_pw.items():
        if mr.upper() == r['receptor'].upper() and ml.upper() == r['ligand'].upper():
            pw = str(p)
            break
    table_data.append([
        i + 1,
        f"{r['receptor']}/{r['ligand']}",
        f"{r['priority_score']:.3f}",
        cons_label if cons_label else '-',
        f"{r['druggability']:.2f}",
        pw if pw else '-',
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
    if 'FGFR2' in str(r['receptor']).upper():
        for j in range(6):
            table[i + 1, j].set_facecolor('#FDEBD0')

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'fig5_cross_species.png', dpi=300, bbox_inches='tight')
plt.savefig(OUTPUT_DIR / 'fig5_cross_species.pdf', bbox_inches='tight')
print("Saved: fig5_cross_species.png/pdf")
