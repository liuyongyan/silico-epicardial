#!/usr/bin/env python3
"""
Figure 1: Cell State Landscape and FGF Family Expression (Mouse)

Panel A: UMAP - Cell States (Normal vs MI)
Panel B: UMAP - FGFR2 Expression
Panel C: UMAP - FGF10 Expression
Panel D: Violin Plots - FGF Family by Cell State
"""

import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import numpy as np
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent.parent.parent
OUTPUT_DIR = PROJECT_DIR / "results" / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================================
# Load mouse data (24GB — use backed mode, then load needed parts)
# ============================================================================
print("Loading mouse data (backed mode)...")
h5ad_path = PROJECT_DIR / "data" / "processed" / "mouse_quaife_ryan_analyzed.h5ad"
adata = sc.read_h5ad(h5ad_path, backed='r')

print(f"Cells: {adata.n_obs}, Genes: {adata.n_vars}")
print(f"obs columns: {list(adata.obs.columns)}")

# Check what metadata is available
print(f"Unique cell states: {adata.obs['cell_state'].unique() if 'cell_state' in adata.obs else 'N/A'}")
if 'condition' in adata.obs:
    print(f"Conditions: {adata.obs['condition'].unique()}")
elif 'disease' in adata.obs:
    print(f"Disease: {adata.obs['disease'].unique()}")

# Check if UMAP exists
has_umap = 'X_umap' in adata.obsm
print(f"Has UMAP: {has_umap}")

if not has_umap:
    print("No UMAP coordinates found. Checking obsm keys...")
    print(f"obsm keys: {list(adata.obsm.keys())}")

# Get UMAP coordinates and metadata into memory
print("Extracting UMAP coordinates and metadata...")
obs_df = adata.obs.copy()

if has_umap:
    umap_coords = adata.obsm['X_umap'].copy() if hasattr(adata.obsm['X_umap'], 'copy') else np.array(adata.obsm['X_umap'])
    obs_df['umap_1'] = umap_coords[:, 0]
    obs_df['umap_2'] = umap_coords[:, 1]

# Get FGF family expression
fgf_genes = ['Fgfr2', 'Fgf10', 'Fgfr1', 'Fgf1', 'Fgf2', 'Fgf7']
available_genes = [g for g in fgf_genes if g in adata.var_names]
print(f"Available FGF genes: {available_genes}")

if available_genes:
    print("Extracting FGF expression...")
    # For backed mode, we need to get expression for specific genes
    gene_indices = sorted([list(adata.var_names).index(g) for g in available_genes])
    # Reorder available_genes to match sorted indices
    idx_to_gene = {list(adata.var_names).index(g): g for g in available_genes}
    available_genes = [idx_to_gene[i] for i in gene_indices]
    # Read the expression data for these genes
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

# Determine condition column
cond_col = None
for col in ['condition', 'disease', 'sample_type']:
    if col in obs_df.columns:
        cond_col = col
        break

if cond_col is None:
    # Try to infer from cell_state
    print("No condition column found, using cell_state only")

# ============================================================================
# Plot
# ============================================================================
if has_umap:
    fig, axes = plt.subplots(2, 2, figsize=(14, 14))

    # Subsample for faster plotting if needed
    n_cells = len(obs_df)
    if n_cells > 50000:
        sample_idx = np.random.RandomState(42).choice(n_cells, 50000, replace=False)
        plot_df = obs_df.iloc[sample_idx].copy()
    else:
        plot_df = obs_df.copy()

    # ---- Panel A: UMAP Cell States ----
    ax = axes[0, 0]
    state_colors = {'quiescent': '#3498DB', 'activated': '#E74C3C'}
    if 'cell_state' in plot_df.columns:
        for state, color in state_colors.items():
            mask = plot_df['cell_state'] == state
            ax.scatter(plot_df.loc[mask, 'umap_1'], plot_df.loc[mask, 'umap_2'],
                       c=color, s=1, alpha=0.3, label=f'{state} (n={mask.sum()})', rasterized=True)
        ax.legend(fontsize=8, markerscale=5)
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    ax.set_title('A. Cell States (Mouse Epicardial)', fontsize=12, fontweight='bold')

    # ---- Panel B: UMAP FGFR2 ----
    ax = axes[0, 1]
    if 'Fgfr2' in plot_df.columns:
        sc_plot = ax.scatter(plot_df['umap_1'], plot_df['umap_2'],
                             c=plot_df['Fgfr2'], cmap='Reds', s=1, alpha=0.5,
                             vmin=0, vmax=np.percentile(plot_df['Fgfr2'][plot_df['Fgfr2'] > 0], 95) if (plot_df['Fgfr2'] > 0).any() else 1,
                             rasterized=True)
        plt.colorbar(sc_plot, ax=ax, shrink=0.6, label='Expression')
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    ax.set_title('B. FGFR2 Expression', fontsize=12, fontweight='bold')

    # ---- Panel C: UMAP FGF10 ----
    ax = axes[1, 0]
    if 'Fgf10' in plot_df.columns:
        sc_plot = ax.scatter(plot_df['umap_1'], plot_df['umap_2'],
                             c=plot_df['Fgf10'], cmap='Blues', s=1, alpha=0.5,
                             vmin=0, vmax=np.percentile(plot_df['Fgf10'][plot_df['Fgf10'] > 0], 95) if (plot_df['Fgf10'] > 0).any() else 1,
                             rasterized=True)
        plt.colorbar(sc_plot, ax=ax, shrink=0.6, label='Expression')
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    ax.set_title('C. FGF10 Expression', fontsize=12, fontweight='bold')

    # ---- Panel D: Violin Plots FGF Family ----
    ax = axes[1, 1]
    if 'cell_state' in plot_df.columns:
        violin_genes = [g for g in ['Fgf1', 'Fgf2', 'Fgf7', 'Fgf10', 'Fgfr1', 'Fgfr2'] if g in plot_df.columns]

        positions = []
        data_q = []
        data_a = []
        labels = []

        for i, gene in enumerate(violin_genes):
            q_vals = plot_df.loc[plot_df['cell_state'] == 'quiescent', gene].values
            a_vals = plot_df.loc[plot_df['cell_state'] == 'activated', gene].values
            data_q.append(q_vals)
            data_a.append(a_vals)
            labels.append(gene)

        x = np.arange(len(labels))
        width = 0.35

        vp1 = ax.violinplot(data_q, positions=x - width/2, showmeans=True, showmedians=False, widths=0.3)
        vp2 = ax.violinplot(data_a, positions=x + width/2, showmeans=True, showmedians=False, widths=0.3)

        for pc in vp1['bodies']:
            pc.set_facecolor('#3498DB')
            pc.set_alpha(0.6)
        for pc in vp2['bodies']:
            pc.set_facecolor('#E74C3C')
            pc.set_alpha(0.6)

        # Style lines
        for partname in ['cmeans', 'cmins', 'cmaxes', 'cbars']:
            if partname in vp1:
                vp1[partname].set_color('#3498DB')
            if partname in vp2:
                vp2[partname].set_color('#E74C3C')

        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=9)
        ax.set_ylabel('Expression (log1p)')
        ax.set_title('D. FGF Family by Cell State', fontsize=12, fontweight='bold')

        # Legend
        from matplotlib.patches import Patch
        ax.legend([Patch(color='#3498DB', alpha=0.6), Patch(color='#E74C3C', alpha=0.6)],
                  ['Quiescent', 'Activated'], fontsize=9)

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'fig1_cell_states_fgf.png', dpi=200, bbox_inches='tight')
    plt.savefig(OUTPUT_DIR / 'fig1_cell_states_fgf.pdf', bbox_inches='tight')
    print("Saved: fig1_cell_states_fgf.png/pdf")
else:
    print("ERROR: No UMAP coordinates in h5ad file. Cannot generate Figure 1.")
    print("Available obsm keys:", list(adata.obsm.keys()) if hasattr(adata, 'obsm') else 'none')
