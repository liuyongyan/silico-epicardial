#!/usr/bin/env python3
"""
Step 5: Compare Results Across Methods

Compares receptor rankings from:
1. Round 1: Embedding perturbation (pre-trained)
2. Round 2: Classifier perturbation (fine-tuned)
3. Original DEG analysis (logFC and padj)

Identifies consistently important receptors and highlights FGFR2.
"""

import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr

# Paths
PROJECT_DIR = Path(__file__).parent.parent.parent
RESULTS_DIR = PROJECT_DIR / "results/geneformer"
DEG_DIR = PROJECT_DIR / "results/deg"


def load_results():
    """
    Load results from all analyses.
    """
    results = {}

    # Round 1
    r1_file = RESULTS_DIR / "perturbation_results_round1.csv"
    if r1_file.exists():
        results['round1'] = pd.read_csv(r1_file)
        print(f"Loaded Round 1: {len(results['round1'])} receptors")
    else:
        print("Round 1 results not found")

    # Round 2
    r2_file = RESULTS_DIR / "perturbation_results_round2.csv"
    if r2_file.exists():
        results['round2'] = pd.read_csv(r2_file)
        print(f"Loaded Round 2: {len(results['round2'])} receptors")
    else:
        print("Round 2 results not found")

    # Original DEG
    deg_file = DEG_DIR / "receptor_deg_LIANA.csv"
    if deg_file.exists():
        deg = pd.read_csv(deg_file)
        deg = deg[(deg['pvals_adj'] < 0.05) & (deg['logfoldchanges'] > 0)]
        deg = deg.rename(columns={'gene': 'gene_symbol', 'names': 'ensembl_id',
                                   'logfoldchanges': 'logFC', 'pvals_adj': 'padj'})
        results['deg'] = deg
        print(f"Loaded DEG: {len(results['deg'])} receptors")

    return results


def create_combined_ranking(results):
    """
    Create a combined ranking across all methods.
    """
    # Get all gene symbols
    all_genes = set()
    for key, df in results.items():
        if 'gene_symbol' in df.columns:
            all_genes.update(df['gene_symbol'].tolist())

    print(f"\nTotal receptors: {len(all_genes)}")

    # Create combined dataframe
    combined = pd.DataFrame({'gene_symbol': list(all_genes)})

    # Add DEG metrics
    if 'deg' in results:
        deg = results['deg'][['gene_symbol', 'logFC', 'padj']].copy()
        deg = deg.rename(columns={'logFC': 'deg_logFC', 'padj': 'deg_padj'})
        combined = combined.merge(deg, on='gene_symbol', how='left')

        # Rank by logFC (descending)
        combined['rank_logFC'] = combined['deg_logFC'].rank(ascending=False)
        # Rank by padj (ascending - lower is better)
        combined['rank_padj'] = combined['deg_padj'].rank(ascending=True)

    # Add Round 1 metrics
    if 'round1' in results:
        r1 = results['round1'].copy()

        # Determine the main metric column
        if 'mean_effect' in r1.columns:
            r1 = r1[['gene_symbol', 'mean_effect']].rename(
                columns={'mean_effect': 'r1_effect'})
        elif 'expression_ratio' in r1.columns:
            r1 = r1[['gene_symbol', 'expression_ratio']].rename(
                columns={'expression_ratio': 'r1_effect'})
        else:
            r1 = r1[['gene_symbol']].copy()
            r1['r1_effect'] = 0

        combined = combined.merge(r1, on='gene_symbol', how='left')
        combined['rank_r1'] = combined['r1_effect'].rank(ascending=False)

    # Add Round 2 metrics
    if 'round2' in results:
        r2 = results['round2'][['gene_symbol', 'mean_prob_drop', 'flip_rate']].copy()
        r2 = r2.rename(columns={'mean_prob_drop': 'r2_prob_drop', 'flip_rate': 'r2_flip_rate'})
        combined = combined.merge(r2, on='gene_symbol', how='left')
        combined['rank_r2'] = combined['r2_prob_drop'].rank(ascending=False)

    # Calculate average rank
    rank_cols = [c for c in combined.columns if c.startswith('rank_')]
    if rank_cols:
        combined['avg_rank'] = combined[rank_cols].mean(axis=1)
        combined['combined_rank'] = combined['avg_rank'].rank()

    return combined


def analyze_fgfr2(combined):
    """
    Detailed analysis of FGFR2's ranking.
    """
    print("\n" + "=" * 60)
    print("FGFR2 Analysis")
    print("=" * 60)

    fgfr2 = combined[combined['gene_symbol'] == 'FGFR2']

    if len(fgfr2) == 0:
        print("FGFR2 not found in results")
        return

    fgfr2 = fgfr2.iloc[0]

    print(f"\nFGFR2 Rankings:")
    print(f"  DEG logFC rank: {fgfr2.get('rank_logFC', 'N/A'):.0f}/84")
    print(f"  DEG padj rank: {fgfr2.get('rank_padj', 'N/A'):.0f}/84")

    if 'rank_r1' in fgfr2.index:
        print(f"  Round 1 (embedding) rank: {fgfr2['rank_r1']:.0f}/84")

    if 'rank_r2' in fgfr2.index:
        print(f"  Round 2 (classifier) rank: {fgfr2['rank_r2']:.0f}/84")

    if 'combined_rank' in fgfr2.index:
        print(f"\n  Combined rank: {fgfr2['combined_rank']:.0f}/84")

    print(f"\nFGFR2 Metrics:")
    print(f"  logFC: {fgfr2.get('deg_logFC', 'N/A'):.3f}")
    print(f"  padj: {fgfr2.get('deg_padj', 'N/A'):.2e}")

    if 'r1_effect' in fgfr2.index:
        print(f"  Round 1 effect: {fgfr2['r1_effect']:.4f}")

    if 'r2_prob_drop' in fgfr2.index:
        print(f"  Round 2 prob drop: {fgfr2['r2_prob_drop']:.4f}")
        print(f"  Round 2 flip rate: {fgfr2['r2_flip_rate']:.4f}")


def create_visualizations(combined, results):
    """
    Create comparison visualizations.
    """
    print("\n" + "-" * 60)
    print("Creating visualizations...")

    fig = plt.figure(figsize=(16, 12))

    # Plot 1: Combined ranking - Top 20
    ax1 = fig.add_subplot(2, 2, 1)
    top20 = combined.nsmallest(20, 'combined_rank') if 'combined_rank' in combined.columns else combined.head(20)
    colors = ['red' if g == 'FGFR2' else ('orange' if g == 'FGFR1' else 'steelblue')
              for g in top20['gene_symbol']]

    y_pos = range(len(top20))
    ax1.barh(y_pos, top20['combined_rank'].values if 'combined_rank' in top20.columns else range(len(top20)),
             color=colors)
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(top20['gene_symbol'])
    ax1.invert_yaxis()
    ax1.set_xlabel('Combined Rank (lower = more important)')
    ax1.set_title('Top 20 Receptors by Combined Ranking')

    # Plot 2: Rank comparison heatmap
    ax2 = fig.add_subplot(2, 2, 2)
    rank_cols = [c for c in combined.columns if c.startswith('rank_') and c != 'combined_rank']

    if len(rank_cols) > 1:
        # Get top 30 by combined rank for readability
        top30 = combined.nsmallest(30, 'combined_rank') if 'combined_rank' in combined.columns else combined.head(30)
        rank_matrix = top30.set_index('gene_symbol')[rank_cols]

        # Rename columns for clarity
        rename_map = {
            'rank_logFC': 'DEG\n(logFC)',
            'rank_padj': 'DEG\n(padj)',
            'rank_r1': 'Round 1\n(Embedding)',
            'rank_r2': 'Round 2\n(Classifier)'
        }
        rank_matrix = rank_matrix.rename(columns=rename_map)

        sns.heatmap(rank_matrix, cmap='RdYlGn_r', annot=True, fmt='.0f',
                    ax=ax2, cbar_kws={'label': 'Rank'})
        ax2.set_title('Ranking Comparison (Top 30)')

        # Highlight FGFR2 row
        fgfr2_idx = list(rank_matrix.index).index('FGFR2') if 'FGFR2' in rank_matrix.index else None
        if fgfr2_idx is not None:
            ax2.add_patch(plt.Rectangle((0, fgfr2_idx), len(rank_cols), 1,
                                         fill=False, edgecolor='red', linewidth=3))

    # Plot 3: Round 1 vs Round 2 comparison
    ax3 = fig.add_subplot(2, 2, 3)
    if 'rank_r1' in combined.columns and 'rank_r2' in combined.columns:
        valid = combined.dropna(subset=['rank_r1', 'rank_r2'])
        ax3.scatter(valid['rank_r1'], valid['rank_r2'], alpha=0.6, s=50, c='steelblue')

        # Highlight FGFR2 and FGFR1
        fgfr2 = combined[combined['gene_symbol'] == 'FGFR2']
        fgfr1 = combined[combined['gene_symbol'] == 'FGFR1']

        if len(fgfr2) > 0:
            ax3.scatter(fgfr2['rank_r1'], fgfr2['rank_r2'],
                       color='red', s=150, zorder=5, label='FGFR2', marker='*')
        if len(fgfr1) > 0:
            ax3.scatter(fgfr1['rank_r1'], fgfr1['rank_r2'],
                       color='orange', s=150, zorder=5, label='FGFR1', marker='*')

        # Add diagonal line
        max_rank = max(valid['rank_r1'].max(), valid['rank_r2'].max())
        ax3.plot([1, max_rank], [1, max_rank], 'k--', alpha=0.3)

        ax3.set_xlabel('Round 1 Rank (Embedding)')
        ax3.set_ylabel('Round 2 Rank (Classifier)')
        ax3.set_title('Round 1 vs Round 2 Ranking')
        ax3.legend()

        # Calculate correlation
        corr, pval = spearmanr(valid['rank_r1'], valid['rank_r2'])
        ax3.text(0.05, 0.95, f'Spearman r = {corr:.3f}\np = {pval:.2e}',
                transform=ax3.transAxes, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Plot 4: DEG rank vs Perturbation rank
    ax4 = fig.add_subplot(2, 2, 4)
    if 'rank_logFC' in combined.columns and 'rank_r2' in combined.columns:
        valid = combined.dropna(subset=['rank_logFC', 'rank_r2'])
        ax4.scatter(valid['rank_logFC'], valid['rank_r2'], alpha=0.6, s=50, c='steelblue')

        # Highlight FGFR2 and FGFR1
        if len(fgfr2) > 0:
            ax4.scatter(fgfr2['rank_logFC'], fgfr2['rank_r2'],
                       color='red', s=150, zorder=5, label='FGFR2', marker='*')
        if len(fgfr1) > 0:
            ax4.scatter(fgfr1['rank_logFC'], fgfr1['rank_r2'],
                       color='orange', s=150, zorder=5, label='FGFR1', marker='*')

        ax4.set_xlabel('DEG Rank (logFC)')
        ax4.set_ylabel('Round 2 Rank (Classifier)')
        ax4.set_title('Expression Change vs Functional Importance')
        ax4.legend()

        # Add annotation
        ax4.text(0.95, 0.05, 'Genes in upper-left:\nLow expression change\nbut high functional importance',
                transform=ax4.transAxes, verticalalignment='bottom', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.3), fontsize=9)

    plt.tight_layout()
    fig.savefig(RESULTS_DIR / "combined_analysis.png", dpi=150, bbox_inches='tight')
    print(f"Saved figure to {RESULTS_DIR / 'combined_analysis.png'}")


def main():
    print("=" * 60)
    print("Comparing Results Across All Methods")
    print("=" * 60)

    # Load results
    results = load_results()

    if len(results) == 0:
        print("No results found. Please run the analysis scripts first.")
        return

    # Create combined ranking
    combined = create_combined_ranking(results)

    # Sort by combined rank
    if 'combined_rank' in combined.columns:
        combined = combined.sort_values('combined_rank')

    # Save combined results
    combined.to_csv(RESULTS_DIR / "combined_rankings.csv", index=False)
    print(f"\nSaved combined rankings to {RESULTS_DIR / 'combined_rankings.csv'}")

    # Display top 20
    print("\n" + "=" * 60)
    print("Top 20 Receptors by Combined Ranking")
    print("=" * 60)

    display_cols = ['gene_symbol', 'combined_rank', 'rank_logFC', 'rank_r1', 'rank_r2', 'deg_logFC']
    display_cols = [c for c in display_cols if c in combined.columns]
    print(combined[display_cols].head(20).to_string(index=False))

    # FGFR2 analysis
    analyze_fgfr2(combined)

    # Compare FGFR1 and FGFR2
    print("\n" + "-" * 60)
    print("FGFR1 vs FGFR2 Comparison")
    print("-" * 60)

    for gene in ['FGFR1', 'FGFR2']:
        row = combined[combined['gene_symbol'] == gene]
        if len(row) > 0:
            row = row.iloc[0]
            print(f"\n{gene}:")
            for col in display_cols:
                if col in row.index and col != 'gene_symbol':
                    val = row[col]
                    if isinstance(val, float):
                        print(f"  {col}: {val:.2f}")
                    else:
                        print(f"  {col}: {val}")

    # Create visualizations
    create_visualizations(combined, results)

    # Final summary
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)

    if 'combined_rank' in combined.columns:
        fgfr2_rank = combined[combined['gene_symbol'] == 'FGFR2']['combined_rank'].values
        if len(fgfr2_rank) > 0:
            print(f"\nFGFR2 Combined Rank: {fgfr2_rank[0]:.0f}/84")

            if fgfr2_rank[0] <= 20:
                print("\n✓ FGFR2 is in the TOP 20 most important receptors!")
            elif fgfr2_rank[0] <= 40:
                print("\n→ FGFR2 is in the TOP HALF of important receptors")
            else:
                print("\n→ FGFR2 ranked in lower half - may need alternative analysis")

    print("\n" + "=" * 60)
    print("Analysis Complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
