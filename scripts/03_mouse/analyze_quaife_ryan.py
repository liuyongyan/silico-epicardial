#!/usr/bin/env python3
"""
Analyze Quaife-Ryan 2021 mouse epicardial dataset
Goal: Classify activation states and analyze FGF10/FGFR2 expression
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.mixture import GaussianMixture
import os

# Paths
OUTPUT_DIR = "/Users/yliu/Desktop/Columbia - Biostatistics/Cheng Lab/Silico Epicardial/data/processed"
RESULTS_DIR = "/Users/yliu/Desktop/Columbia - Biostatistics/Cheng Lab/Silico Epicardial/results/mouse"

os.makedirs(RESULTS_DIR, exist_ok=True)
sc.settings.figdir = RESULTS_DIR

print("=" * 60)
print("Analyzing Quaife-Ryan 2021 Mouse Epicardial Dataset")
print("=" * 60)

# Load data
print("\n1. Loading processed data...")
adata = sc.read(f"{OUTPUT_DIR}/mouse_quaife_ryan_raw.h5ad")
print(f"   Shape: {adata.shape}")
print(f"   Disease distribution:\n{adata.obs['disease'].value_counts()}")

# This dataset is specifically epicardial stromal cells - all cells are epicardial-derived
# According to the paper, they used a novel perfusion method to isolate epicardial cells
print("\n2. Dataset note: All cells are epicardial stromal cells (EpiSC) by design")
print("   - Paper used perfusion-based isolation specifically for epicardial cells")
print("   - No need to filter for epicardial markers")

# Normalize if not already (data appears to be TPM-like normalized)
print("\n3. Preparing data for analysis...")
# Store raw counts
adata.raw = adata.copy()

# Log transform (data is normalized but not log-transformed)
adata.X = np.log1p(adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X)
print("   Log-transformed expression values")

# Score cells for activation state
print("\n4. Scoring cells for activation state...")

# Quiescent/epithelial signature
quiescent_genes = ['Wt1', 'Upk3b', 'Msln', 'Krt19', 'Cdh1', 'Cldn1', 'Cldn3', 'Tjp1', 'Ocln']
quiescent_in_data = [g for g in quiescent_genes if g in adata.var_names]
print(f"   Quiescent genes used: {quiescent_in_data}")

# Activated/mesenchymal signature
activated_genes = ['Snai1', 'Snai2', 'Twist1', 'Vim', 'Fn1', 'Col1a1', 'Col3a1', 'Postn', 'Acta2']
activated_in_data = [g for g in activated_genes if g in adata.var_names]
print(f"   Activated genes used: {activated_in_data}")

sc.tl.score_genes(adata, gene_list=quiescent_in_data, score_name='quiescent_score', ctrl_size=50)
sc.tl.score_genes(adata, gene_list=activated_in_data, score_name='activated_score', ctrl_size=50)

# Calculate state potential
adata.obs['state_potential'] = adata.obs['activated_score'] - adata.obs['quiescent_score']

# Use GMM to classify cells
print("\n5. Classifying cells using GMM...")
state_values = adata.obs['state_potential'].values.reshape(-1, 1)
gmm = GaussianMixture(n_components=2, random_state=42)
gmm.fit(state_values)
labels = gmm.predict(state_values)

# Determine which label is activated (higher mean)
means = gmm.means_.flatten()
activated_label = np.argmax(means)
quiescent_label = 1 - activated_label

adata.obs['cell_state'] = ['activated' if l == activated_label else 'quiescent' for l in labels]
print(f"   Cell state distribution:\n{adata.obs['cell_state'].value_counts()}")

# Also identify transitioning cells based on probability
probs = gmm.predict_proba(state_values)
adata.obs['activation_prob'] = probs[:, activated_label]
adata.obs['cell_state_detailed'] = adata.obs['cell_state'].copy()
transitioning = (adata.obs['activation_prob'] > 0.3) & (adata.obs['activation_prob'] < 0.7)
adata.obs.loc[transitioning, 'cell_state_detailed'] = 'transitioning'
print(f"   Detailed state distribution:\n{adata.obs['cell_state_detailed'].value_counts()}")

# Analyze FGF expression
print("\n6. Analyzing FGF family expression...")
fgf_genes = ['Fgf1', 'Fgf2', 'Fgf7', 'Fgf10', 'Fgfr1', 'Fgfr2', 'Fgfr3']
fgf_in_data = [g for g in fgf_genes if g in adata.var_names]

# Get expression by cell state
fgf_expr = pd.DataFrame()
for gene in fgf_in_data:
    gene_idx = adata.var_names.tolist().index(gene)
    expr = adata.X[:, gene_idx] if not hasattr(adata.X, 'toarray') else adata.X[:, gene_idx].toarray().flatten()
    fgf_expr[gene] = expr

fgf_expr['cell_state'] = adata.obs['cell_state'].values

# Calculate mean expression by state
print("\n   FGF expression by cell state:")
print("   " + "-" * 50)
print(f"   {'Gene':<10} {'Quiescent':>12} {'Activated':>12} {'Log2FC':>10}")
print("   " + "-" * 50)

results = []
for gene in fgf_in_data:
    quiescent_expr = fgf_expr[fgf_expr['cell_state'] == 'quiescent'][gene].mean()
    activated_expr = fgf_expr[fgf_expr['cell_state'] == 'activated'][gene].mean()

    # Avoid log(0) by adding pseudocount
    log2fc = np.log2((activated_expr + 0.01) / (quiescent_expr + 0.01))

    print(f"   {gene:<10} {quiescent_expr:>12.3f} {activated_expr:>12.3f} {log2fc:>10.2f}")
    results.append({
        'gene': gene,
        'quiescent_mean': quiescent_expr,
        'activated_mean': activated_expr,
        'log2fc': log2fc
    })

# Differential expression analysis
print("\n7. Running differential expression analysis...")
sc.tl.rank_genes_groups(adata, groupby='cell_state', method='wilcoxon', key_added='de_state')

# Get DE results for activated vs quiescent
de_results = sc.get.rank_genes_groups_df(adata, group='activated', key='de_state')
print(f"   Top upregulated in activated cells:")
print(de_results.head(20).to_string())

# Check FGF genes in DE results
print("\n   FGF genes in DE results:")
for gene in fgf_in_data:
    gene_de = de_results[de_results['names'] == gene]
    if len(gene_de) > 0:
        print(f"   {gene}: logFC={gene_de['logfoldchanges'].values[0]:.3f}, pval_adj={gene_de['pvals_adj'].values[0]:.2e}")

# Create visualizations
print("\n8. Creating visualizations...")

# Run PCA and UMAP
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat')
sc.pp.pca(adata, use_highly_variable=True)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# Plot UMAP with cell states
fig, axes = plt.subplots(2, 3, figsize=(15, 10))

# UMAP by cell state
sc.pl.umap(adata, color='cell_state', ax=axes[0, 0], show=False, title='Cell State')
sc.pl.umap(adata, color='activation_prob', ax=axes[0, 1], show=False, title='Activation Probability', cmap='RdYlBu_r')
sc.pl.umap(adata, color='disease', ax=axes[0, 2], show=False, title='Disease Condition')

# UMAP with key genes
sc.pl.umap(adata, color='Fgf10', ax=axes[1, 0], show=False, title='Fgf10', cmap='viridis')
sc.pl.umap(adata, color='Fgfr2', ax=axes[1, 1], show=False, title='Fgfr2', cmap='viridis')
sc.pl.umap(adata, color='Wt1', ax=axes[1, 2], show=False, title='Wt1 (epicardial marker)', cmap='viridis')

plt.tight_layout()
plt.savefig(f"{RESULTS_DIR}/quaife_ryan_umap_overview.png", dpi=150, bbox_inches='tight')
plt.close()
print(f"   Saved: {RESULTS_DIR}/quaife_ryan_umap_overview.png")

# FGF expression violin plot
fig, ax = plt.subplots(figsize=(12, 6))
fgf_melted = fgf_expr.melt(id_vars=['cell_state'], value_vars=fgf_in_data)
sns.violinplot(data=fgf_melted, x='variable', y='value', hue='cell_state', split=True, ax=ax)
ax.set_xlabel('Gene')
ax.set_ylabel('Log Expression')
ax.set_title('FGF Family Expression by Cell State (Quaife-Ryan 2021)')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(f"{RESULTS_DIR}/quaife_ryan_fgf_violin.png", dpi=150, bbox_inches='tight')
plt.close()
print(f"   Saved: {RESULTS_DIR}/quaife_ryan_fgf_violin.png")

# Save results
print("\n9. Saving results...")
adata.write(f"{OUTPUT_DIR}/mouse_quaife_ryan_analyzed.h5ad")
print(f"   Saved: {OUTPUT_DIR}/mouse_quaife_ryan_analyzed.h5ad")

# Save DE results
de_results.to_csv(f"{RESULTS_DIR}/quaife_ryan_de_activated_vs_quiescent.csv", index=False)
print(f"   Saved: {RESULTS_DIR}/quaife_ryan_de_activated_vs_quiescent.csv")

# Save FGF summary
fgf_summary = pd.DataFrame(results)
fgf_summary.to_csv(f"{RESULTS_DIR}/quaife_ryan_fgf_summary.csv", index=False)
print(f"   Saved: {RESULTS_DIR}/quaife_ryan_fgf_summary.csv")

print("\n" + "=" * 60)
print("Analysis complete!")
print("=" * 60)
print(f"\nKey findings:")
print(f"  - Quiescent cells: {(adata.obs['cell_state'] == 'quiescent').sum():,}")
print(f"  - Activated cells: {(adata.obs['cell_state'] == 'activated').sum():,}")
print(f"  - Transitioning: {(adata.obs['cell_state_detailed'] == 'transitioning').sum():,}")
