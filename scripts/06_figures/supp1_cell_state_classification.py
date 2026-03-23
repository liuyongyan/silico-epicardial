#!/usr/bin/env python3
"""
Supplementary Figure 1: Cell State Classification Method (Mouse)

Panel A: Signature Gene Expression Heatmap (quiescent vs activated markers)
Panel B: GMM Threshold Selection (state_potential histogram + GMM fit)
Panel C: Activation Probability Distribution by Condition (Normal vs MI)
"""

import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent.parent.parent
OUTPUT_DIR = PROJECT_DIR / "results" / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================================
# Load data
# ============================================================================
print("Loading mouse data (backed mode)...")
adata = sc.read_h5ad(PROJECT_DIR / "data" / "processed" / "mouse_quaife_ryan_analyzed.h5ad", backed='r')
obs_df = adata.obs.copy()
print(f"Cells: {len(obs_df)}")

# Signature genes
quiescent_genes = ['Wt1', 'Upk3b', 'Msln', 'Krt19', 'Cdh1', 'Cldn1', 'Cldn3', 'Tjp1', 'Ocln']
activated_genes = ['Snai1', 'Snai2', 'Twist1', 'Vim', 'Fn1', 'Col1a1', 'Col3a1', 'Postn', 'Acta2']
all_sig_genes = quiescent_genes + activated_genes

# Get expression for signature genes
available = [g for g in all_sig_genes if g in adata.var_names]
print(f"Signature genes available: {len(available)}/{len(all_sig_genes)}")

gene_indices = sorted([list(adata.var_names).index(g) for g in available])
idx_to_gene = {list(adata.var_names).index(g): g for g in available}
available_sorted = [idx_to_gene[i] for i in gene_indices]

print("Extracting signature gene expression...")
expr = adata.X[:, gene_indices]
if hasattr(expr, 'toarray'):
    expr = expr.toarray()
elif hasattr(expr, 'todense'):
    expr = np.asarray(expr.todense())
else:
    expr = np.array(expr)

expr_df = pd.DataFrame(expr, columns=available_sorted, index=obs_df.index)
print("Data loaded.")

# ============================================================================
# Plot
# ============================================================================
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# ---- Panel A: Signature Gene Heatmap ----
ax = axes[0]

# Compute mean expression per cell state for each signature gene
mean_expr = {}
for state in ['quiescent', 'activated']:
    mask = obs_df['cell_state'] == state
    mean_expr[state] = expr_df.loc[mask].mean()

hm_df = pd.DataFrame(mean_expr)
# Reorder: quiescent genes first, then activated genes
gene_order = [g for g in quiescent_genes if g in hm_df.index] + \
             [g for g in activated_genes if g in hm_df.index]
hm_df = hm_df.loc[gene_order]

import seaborn as sns

# Z-score normalize per gene (row) for better visualization
hm_z = hm_df.subtract(hm_df.mean(axis=1), axis=0).divide(hm_df.std(axis=1) + 1e-10, axis=0)

sns.heatmap(hm_z, cmap='RdBu_r', center=0, annot=hm_df.round(2), fmt='',
            annot_kws={'fontsize': 8}, ax=ax,
            cbar_kws={'shrink': 0.6, 'label': 'Z-score'},
            yticklabels=True)

# Add divider between quiescent and activated genes
n_q = len([g for g in quiescent_genes if g in hm_df.index])
ax.axhline(n_q, color='black', linewidth=2)
ax.text(-0.5, n_q/2, 'Quiescent\nmarkers', fontsize=8, ha='right', va='center',
        fontweight='bold', color='#3498DB', rotation=0)
ax.text(-0.5, n_q + (len(gene_order)-n_q)/2, 'Activated\nmarkers', fontsize=8,
        ha='right', va='center', fontweight='bold', color='#E74C3C', rotation=0)

ax.set_title('A. Signature Gene Expression', fontsize=12, fontweight='bold')

# ---- Panel B: GMM Threshold on state_potential ----
ax = axes[1]

sp = obs_df['state_potential'].values
# Clip extremes for visualization
sp_clip = np.clip(sp, np.percentile(sp, 1), np.percentile(sp, 99))

# Fit GMM
gmm = GaussianMixture(n_components=2, random_state=42)
gmm.fit(sp_clip.reshape(-1, 1))
means = gmm.means_.flatten()
stds = np.sqrt(gmm.covariances_.flatten())
weights = gmm.weights_

threshold = np.mean(means)

# Histogram
ax.hist(sp_clip, bins=100, density=True, alpha=0.5, color='#95A5A6', edgecolor='none')

# GMM components
x_range = np.linspace(sp_clip.min(), sp_clip.max(), 500)
from scipy.stats import norm
for i, (mean, std, w, color, label) in enumerate(zip(
    means, stds, weights,
    ['#3498DB', '#E74C3C'],
    ['Low (quiescent)', 'High (activated)']
)):
    ax.plot(x_range, w * norm.pdf(x_range, mean, std), color=color, linewidth=2, label=label)

# Threshold line
ax.axvline(threshold, color='black', linewidth=2, linestyle='--', label=f'Threshold = {threshold:.0f}')

ax.text(threshold + 0.5, ax.get_ylim()[1] * 0.5 if ax.get_ylim()[1] > 0 else 0.5,
        'Agreement: 60.8%', fontsize=9, color='gray', style='italic',
        ha='left', va='center')

ax.set_xlabel('State Potential (activated_score − quiescent_score)', fontsize=10)
ax.set_ylabel('Density', fontsize=10)
ax.set_title('B. GMM Threshold Selection', fontsize=12, fontweight='bold')
ax.legend(fontsize=8)

# ---- Panel C: Activation Probability by Condition ----
ax = axes[2]

for condition, color, label in [
    ('normal', '#3498DB', 'Normal'),
    ('myocardial infarction', '#E74C3C', 'MI')
]:
    mask = obs_df['disease'] == condition
    vals = obs_df.loc[mask, 'activation_prob']
    ax.hist(vals, bins=50, density=True, alpha=0.5, color=color, label=f'{label} (n={mask.sum()})',
            edgecolor='none')

    # Add mean line
    mean_val = vals.mean()
    ax.axvline(mean_val, color=color, linewidth=2, linestyle='--', alpha=0.8)
    ax.text(mean_val + 0.02, ax.get_ylim()[1] * 0.8 if condition == 'normal' else ax.get_ylim()[1] * 0.6,
            f'mean={mean_val:.2f}', fontsize=8, color=color, fontweight='bold')

ax.set_xlabel('Activation Probability', fontsize=10)
ax.set_ylabel('Density', fontsize=10)
ax.set_title('C. Activation Probability by Condition', fontsize=12, fontweight='bold')
ax.legend(fontsize=9)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'supp1_cell_state_classification.png', dpi=300, bbox_inches='tight')
plt.savefig(OUTPUT_DIR / 'supp1_cell_state_classification.pdf', bbox_inches='tight')
print("Saved: supp1_cell_state_classification.png/pdf")
