#!/usr/bin/env python3
"""
Figure 1: Cell State Landscape and FGF Family Expression (Mouse)

Panel A: UMAP - Cell States (quiescent vs activated)
Panel B: UMAP - FGFR2 Expression
Panel C: UMAP - FGF10 Expression
Panel D: Dot Plot - FGF Family by Cell State (pct expressing + mean expression)
"""

import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent.parent.parent
OUTPUT_DIR = PROJECT_DIR / "results" / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================================
# Load mouse data
# ============================================================================
print("Loading mouse data (backed mode)...")
h5ad_path = PROJECT_DIR / "data" / "processed" / "mouse_quaife_ryan_analyzed.h5ad"
adata = sc.read_h5ad(h5ad_path, backed='r')

print(f"Cells: {adata.n_obs}, Genes: {adata.n_vars}")

# Extract metadata
obs_df = adata.obs.copy()

# UMAP
umap_coords = np.array(adata.obsm['X_umap'])
obs_df['umap_1'] = umap_coords[:, 0]
obs_df['umap_2'] = umap_coords[:, 1]

# FGF expression
fgf_genes = ['Fgf1', 'Fgf2', 'Fgf7', 'Fgf10', 'Fgfr1', 'Fgfr2']
available_genes = [g for g in fgf_genes if g in adata.var_names]
gene_indices = sorted([list(adata.var_names).index(g) for g in available_genes])
idx_to_gene = {list(adata.var_names).index(g): g for g in available_genes}
available_genes = [idx_to_gene[i] for i in gene_indices]

print(f"Extracting expression for: {available_genes}")
expr_data = adata.X[:, gene_indices]
if hasattr(expr_data, 'toarray'):
    expr_data = expr_data.toarray()
elif hasattr(expr_data, 'todense'):
    expr_data = np.asarray(expr_data.todense())
else:
    expr_data = np.array(expr_data)

for i, gene in enumerate(available_genes):
    obs_df[gene] = expr_data[:, i]

print("Data loaded. Plotting...")

# Subsample for UMAP panels
n_cells = len(obs_df)
if n_cells > 50000:
    sample_idx = np.random.RandomState(42).choice(n_cells, 50000, replace=False)
    plot_df = obs_df.iloc[sample_idx].copy()
else:
    plot_df = obs_df.copy()

# ============================================================================
# Plot
# ============================================================================
fig = plt.figure(figsize=(16, 14))

# Use gridspec for flexible layout
gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)

# ---- Panel A: UMAP Cell States ----
ax = fig.add_subplot(gs[0, 0])
for state, color, zorder in [('quiescent', '#3498DB', 1), ('activated', '#E74C3C', 2)]:
    mask = plot_df['cell_state'] == state
    ax.scatter(plot_df.loc[mask, 'umap_1'], plot_df.loc[mask, 'umap_2'],
               c=color, s=0.5, alpha=0.3, label=f'{state} (n={mask.sum()})',
               rasterized=True, zorder=zorder)
ax.legend(fontsize=8, markerscale=10, loc='upper left')
ax.set_xlabel('UMAP 1', fontsize=10)
ax.set_ylabel('UMAP 2', fontsize=10)
ax.set_title('A. Cell States (Mouse Epicardial)', fontsize=12, fontweight='bold')

# ---- Panel B: UMAP FGFR2 ----
ax = fig.add_subplot(gs[0, 1])
# Sort by expression so high values plotted on top
order = plot_df['Fgfr2'].argsort()
vmax = np.percentile(plot_df['Fgfr2'][plot_df['Fgfr2'] > 0], 95) if (plot_df['Fgfr2'] > 0).any() else 1
sc_plot = ax.scatter(plot_df.iloc[order]['umap_1'], plot_df.iloc[order]['umap_2'],
                     c=plot_df.iloc[order]['Fgfr2'], cmap='Reds', s=0.5, alpha=0.5,
                     vmin=0, vmax=vmax, rasterized=True)
plt.colorbar(sc_plot, ax=ax, shrink=0.6, label='Expression (log1p)')
ax.set_xlabel('UMAP 1', fontsize=10)
ax.set_ylabel('UMAP 2', fontsize=10)
ax.set_title('B. FGFR2 Expression', fontsize=12, fontweight='bold')

# ---- Panel C: UMAP FGF10 ----
ax = fig.add_subplot(gs[1, 0])
order = plot_df['Fgf10'].argsort()
vmax = np.percentile(plot_df['Fgf10'][plot_df['Fgf10'] > 0], 95) if (plot_df['Fgf10'] > 0).any() else 1
sc_plot = ax.scatter(plot_df.iloc[order]['umap_1'], plot_df.iloc[order]['umap_2'],
                     c=plot_df.iloc[order]['Fgf10'], cmap='Blues', s=0.5, alpha=0.5,
                     vmin=0, vmax=vmax, rasterized=True)
plt.colorbar(sc_plot, ax=ax, shrink=0.6, label='Expression (log1p)')
ax.set_xlabel('UMAP 1', fontsize=10)
ax.set_ylabel('UMAP 2', fontsize=10)
ax.set_title('C. FGF10 Expression', fontsize=12, fontweight='bold')

# ---- Panel D: Violin (all cells incl. zeros, no Fgfr1) + annotations ----
ax = fig.add_subplot(gs[1, 1])

genes_to_plot = ['Fgf1', 'Fgf2', 'Fgf7', 'Fgf10', 'Fgfr2']  # Fgfr1 removed (dominates scale)
genes_to_plot = [g for g in genes_to_plot if g in obs_df.columns]

data_q = []
data_a = []
mean_q = []
mean_a = []
pct_q = []
pct_a = []
for gene in genes_to_plot:
    q_vals = obs_df.loc[obs_df['cell_state'] == 'quiescent', gene]
    a_vals = obs_df.loc[obs_df['cell_state'] == 'activated', gene]
    mean_q.append(q_vals.mean())  # mean over ALL cells (incl. zeros)
    mean_a.append(a_vals.mean())
    pct_q.append((q_vals > 0).mean() * 100)
    pct_a.append((a_vals > 0).mean() * 100)
    data_q.append(q_vals[q_vals > 0].values)  # violin: expressing cells only
    data_a.append(a_vals[a_vals > 0].values)

x = np.arange(len(genes_to_plot))
width = 0.35

vp1 = ax.violinplot(data_q, positions=x - width/2, showmeans=True, showmedians=False, widths=0.3)
vp2 = ax.violinplot(data_a, positions=x + width/2, showmeans=True, showmedians=False, widths=0.3)

for pc in vp1['bodies']:
    pc.set_facecolor('#3498DB')
    pc.set_alpha(0.7)
for pc in vp2['bodies']:
    pc.set_facecolor('#E74C3C')
    pc.set_alpha(0.7)

for partname in ['cmeans', 'cmins', 'cmaxes', 'cbars']:
    if partname in vp1:
        vp1[partname].set_color('#3498DB')
        vp1[partname].set_linewidth(1)
    if partname in vp2:
        vp2[partname].set_color('#E74C3C')
        vp2[partname].set_linewidth(1)

ax.set_xticks(x)
ax.set_xticklabels(genes_to_plot, fontsize=10)
ax.set_ylabel('Expression (log1p)', fontsize=10)
ax.set_title('D. FGF Family by Cell State', fontsize=12, fontweight='bold')

from matplotlib.patches import Patch
ax.legend([Patch(color='#3498DB', alpha=0.7), Patch(color='#E74C3C', alpha=0.7)],
          ['Quiescent', 'Activated'], fontsize=9, loc='upper right')

# Annotations below: % expressing and fold change
for i in range(len(genes_to_plot)):
    fc = mean_a[i] / mean_q[i] if mean_q[i] > 0 else 0
    direction = '↑' if fc > 1.2 else ('↓' if fc < 0.8 else '–')
    fc_color = '#C0392B' if fc > 1.2 else ('#2471A3' if fc < 0.8 else '#7F8C8D')

    # % expressing
    ax.text(x[i] - width/2, -0.6, f'{pct_q[i]:.0f}%', fontsize=8, ha='center', color='#1A5276', fontweight='bold')
    ax.text(x[i] + width/2, -0.6, f'{pct_a[i]:.0f}%', fontsize=8, ha='center', color='#922B21', fontweight='bold')

    # Fold change
    ax.text(x[i], -1.1, f'{fc:.1f}x {direction}', fontsize=10, ha='center', fontweight='bold', color=fc_color)

ax.text(-0.7, -0.6, '% expr:', fontsize=8, ha='right', color='#2C3E50', fontweight='bold')
ax.text(-0.7, -1.1, 'FC:', fontsize=8, ha='right', color='#2C3E50', fontweight='bold')
ax.set_ylim(-1.5, None)

plt.savefig(OUTPUT_DIR / 'fig1_cell_states_fgf.png', dpi=200, bbox_inches='tight')
plt.savefig(OUTPUT_DIR / 'fig1_cell_states_fgf.pdf', bbox_inches='tight')
print("Saved: fig1_cell_states_fgf.png/pdf")
