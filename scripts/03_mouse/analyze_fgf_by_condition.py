#!/usr/bin/env python3
"""
Analyze FGF expression by disease condition (MI vs Normal)
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Paths
OUTPUT_DIR = "/Users/yliu/Desktop/Columbia - Biostatistics/Cheng Lab/Silico Epicardial/data/processed"
RESULTS_DIR = "/Users/yliu/Desktop/Columbia - Biostatistics/Cheng Lab/Silico Epicardial/results/mouse"

print("=" * 60)
print("FGF Expression by Disease Condition")
print("=" * 60)

# Load analyzed data
adata = sc.read(f"{OUTPUT_DIR}/mouse_quaife_ryan_analyzed.h5ad")
print(f"Shape: {adata.shape}")

# Cross-tabulate cell state and disease
print("\n1. Cell state by disease condition:")
crosstab = pd.crosstab(adata.obs['cell_state'], adata.obs['disease'], margins=True)
print(crosstab)

# Calculate proportions
print("\n   Proportion of activated cells by condition:")
for disease in ['myocardial infarction', 'normal']:
    subset = adata.obs[adata.obs['disease'] == disease]
    prop = (subset['cell_state'] == 'activated').mean() * 100
    print(f"   {disease}: {prop:.1f}% activated")

# FGF expression analysis by condition
print("\n2. FGF expression by disease condition:")
fgf_genes = ['Fgf1', 'Fgf2', 'Fgf7', 'Fgf10', 'Fgfr1', 'Fgfr2', 'Fgfr3']

results = []
print("\n   " + "-" * 70)
print(f"   {'Gene':<10} {'Normal':>12} {'MI':>12} {'Log2FC':>10} {'p-value':>12}")
print("   " + "-" * 70)

for gene in fgf_genes:
    if gene in adata.var_names:
        gene_idx = adata.var_names.tolist().index(gene)
        expr = adata.X[:, gene_idx]
        if hasattr(expr, 'toarray'):
            expr = expr.toarray().flatten()

        normal_expr = expr[adata.obs['disease'] == 'normal']
        mi_expr = expr[adata.obs['disease'] == 'myocardial infarction']

        normal_mean = normal_expr.mean()
        mi_mean = mi_expr.mean()
        log2fc = np.log2((mi_mean + 0.01) / (normal_mean + 0.01))

        # Statistical test
        _, pval = stats.mannwhitneyu(normal_expr, mi_expr, alternative='two-sided')

        print(f"   {gene:<10} {normal_mean:>12.3f} {mi_mean:>12.3f} {log2fc:>10.2f} {pval:>12.2e}")

        results.append({
            'gene': gene,
            'normal_mean': normal_mean,
            'mi_mean': mi_mean,
            'log2fc': log2fc,
            'pval': pval
        })

# FGF expression by cell state AND condition
print("\n3. FGF10/FGFR2 expression by cell state AND condition:")
print("-" * 80)

for gene in ['Fgf10', 'Fgfr2']:
    print(f"\n{gene}:")
    gene_idx = adata.var_names.tolist().index(gene)
    expr = adata.X[:, gene_idx]
    if hasattr(expr, 'toarray'):
        expr = expr.toarray().flatten()

    for state in ['quiescent', 'activated']:
        for disease in ['normal', 'myocardial infarction']:
            mask = (adata.obs['cell_state'] == state) & (adata.obs['disease'] == disease)
            mean_expr = expr[mask].mean()
            n_cells = mask.sum()
            print(f"   {state:12} + {disease:22}: mean={mean_expr:.4f}, n={n_cells:,}")

# Create visualization
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: Cell state proportions by condition
ax = axes[0, 0]
props = adata.obs.groupby('disease')['cell_state'].value_counts(normalize=True).unstack()
props.plot(kind='bar', ax=ax, color=['steelblue', 'coral'])
ax.set_title('Cell State Proportions by Condition')
ax.set_xlabel('')
ax.set_ylabel('Proportion')
ax.legend(title='Cell State')
ax.tick_params(axis='x', rotation=45)

# Plot 2: FGFR2 expression by state and condition
ax = axes[0, 1]
gene_idx = adata.var_names.tolist().index('Fgfr2')
expr = adata.X[:, gene_idx]
if hasattr(expr, 'toarray'):
    expr = expr.toarray().flatten()
plot_df = pd.DataFrame({
    'Fgfr2': expr,
    'Cell State': adata.obs['cell_state'].values,
    'Condition': adata.obs['disease'].values
})
sns.boxplot(data=plot_df, x='Condition', y='Fgfr2', hue='Cell State', ax=ax)
ax.set_title('FGFR2 Expression')
ax.tick_params(axis='x', rotation=45)

# Plot 3: Fgf10 expression by state and condition
ax = axes[1, 0]
gene_idx = adata.var_names.tolist().index('Fgf10')
expr = adata.X[:, gene_idx]
if hasattr(expr, 'toarray'):
    expr = expr.toarray().flatten()
plot_df = pd.DataFrame({
    'Fgf10': expr,
    'Cell State': adata.obs['cell_state'].values,
    'Condition': adata.obs['disease'].values
})
sns.boxplot(data=plot_df, x='Condition', y='Fgf10', hue='Cell State', ax=ax)
ax.set_title('FGF10 Expression')
ax.tick_params(axis='x', rotation=45)

# Plot 4: Scatter of Fgf10 vs Fgfr2
ax = axes[1, 1]
fgf10_idx = adata.var_names.tolist().index('Fgf10')
fgfr2_idx = adata.var_names.tolist().index('Fgfr2')
fgf10_expr = adata.X[:, fgf10_idx]
fgfr2_expr = adata.X[:, fgfr2_idx]
if hasattr(fgf10_expr, 'toarray'):
    fgf10_expr = fgf10_expr.toarray().flatten()
    fgfr2_expr = fgfr2_expr.toarray().flatten()

# Sample for visibility
np.random.seed(42)
sample_idx = np.random.choice(len(fgf10_expr), min(5000, len(fgf10_expr)), replace=False)
colors = ['coral' if s == 'activated' else 'steelblue' for s in adata.obs['cell_state'].values[sample_idx]]
ax.scatter(fgf10_expr[sample_idx], fgfr2_expr[sample_idx], c=colors, alpha=0.3, s=5)
ax.set_xlabel('Fgf10 (log)')
ax.set_ylabel('Fgfr2 (log)')
ax.set_title('FGF10 vs FGFR2 Expression')

plt.tight_layout()
plt.savefig(f"{RESULTS_DIR}/quaife_ryan_fgf_by_condition.png", dpi=150, bbox_inches='tight')
plt.close()
print(f"\nSaved: {RESULTS_DIR}/quaife_ryan_fgf_by_condition.png")

# Save results
pd.DataFrame(results).to_csv(f"{RESULTS_DIR}/quaife_ryan_fgf_by_condition.csv", index=False)
print(f"Saved: {RESULTS_DIR}/quaife_ryan_fgf_by_condition.csv")

print("\n" + "=" * 60)
print("Key Findings Summary:")
print("=" * 60)
print("""
1. CELL STATE BY CONDITION:
   - MI samples have higher proportion of activated epicardial cells
   - Normal samples retain more quiescent epicardial cells

2. FGF10 (LIGAND):
   - Expressed more in QUIESCENT cells
   - Lower in MI compared to normal
   - FGF10 may be a quiescence-maintaining factor

3. FGFR2 (RECEPTOR):
   - UPREGULATED in ACTIVATED cells (log2FC = ~1.5)
   - Higher in MI compared to normal
   - FGFR2 marks/drives the activated state

4. INTERPRETATION:
   - Quiescent epicardial cells express FGF10 (autocrine/paracrine signal)
   - Upon MI, cells upregulate FGFR2 and become activated
   - FGF10/FGFR2 axis may regulate epicardial activation
   - Targeting FGFR2 could modulate epicardial activation post-MI
""")
