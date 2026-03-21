#!/usr/bin/env python3
"""
Human L-R Mismatch Analysis (PERIHEART Epicardial Cells)

Same methodology as mouse analysis (02_mouse_lr_mismatch_refined.py):
  composite = r_norm × l_norm × starvation_ratio × db_weight

Uses two DEG comparisons:
  A) Quiescent vs Activated (cell state, matches mouse analysis)
  B) Normal vs MI (condition-based)
"""

import pandas as pd
import numpy as np
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent.parent.parent
DEG_DIR = PROJECT_DIR / "results" / "deg"
OUTPUT_DIR = PROJECT_DIR / "results" / "mismatch"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================================
# 1. Build L-R database with consensus scores (human, no conversion needed)
# ============================================================================
print("=" * 80)
print("Step 1: Loading human L-R database from OmniPath")
print("=" * 80)

import omnipath

lr_raw = omnipath.interactions.PostTranslational.get(genesymbols=True, organism='human')
sources_filter = 'CellPhoneDB|CellTalkDB|Ramilowski2015|Fantom5_LRdb|connectomeDB2020'
lr_curated = lr_raw[lr_raw['sources'].str.contains(sources_filter, na=False)].copy()

target_dbs = ['CellPhoneDB', 'CellTalkDB', 'Ramilowski2015', 'Fantom5_LRdb', 'connectomeDB2020']

def count_dbs(sources_str):
    return sum(1 for db in target_dbs if db in sources_str)

lr_curated['n_db'] = lr_curated['sources'].apply(count_dbs)

# Extract and expand complexes
raw_pairs = lr_curated[['source_genesymbol', 'target_genesymbol', 'n_db']].copy()
raw_pairs.columns = ['ligand', 'receptor', 'n_db']

expanded = []
for _, row in raw_pairs.iterrows():
    ligs = row['ligand'].split('_') if '_' in row['ligand'] else [row['ligand']]
    recs = row['receptor'].split('_') if '_' in row['receptor'] else [row['receptor']]
    for l in ligs:
        for r in recs:
            expanded.append({'ligand': l, 'receptor': r, 'n_db': row['n_db']})

lr_exp = pd.DataFrame(expanded)
lr_exp = lr_exp.groupby(['ligand', 'receptor'])['n_db'].max().reset_index()
print(f"Human L-R pairs: {len(lr_exp)}")

# ============================================================================
# 2. Load human DEG data
# ============================================================================
print("\n" + "=" * 80)
print("Step 2: Loading human DEG data")
print("=" * 80)

# --- A) Full Wilcoxon DEG (quiescent vs activated) ---
deg_wilcox = pd.read_csv(DEG_DIR / "deg_results_wilcoxon.csv")
# Column 'gene' has gene symbol, 'names' has ensembl ID
deg_wilcox = deg_wilcox.dropna(subset=['logfoldchanges'])
deg_wilcox = deg_wilcox[np.isfinite(deg_wilcox['logfoldchanges'])]
print(f"Wilcoxon DEG (all genes): {len(deg_wilcox)}")

# --- B) Receptor-level DEGs ---
rec_state = pd.read_csv(DEG_DIR / "receptor_quiescent_vs_activated.csv")
rec_cond = pd.read_csv(DEG_DIR / "receptor_normal_vs_mi_PERIHEART.csv")
print(f"Receptors (quiescent vs activated): {len(rec_state)}")
print(f"Receptors (normal vs MI): {len(rec_cond)}")

# ============================================================================
# 3. Run mismatch analysis for BOTH comparisons
# ============================================================================

DB_THRESHOLD = 3

def run_mismatch(receptor_df, deg_full_df, comparison_name,
                 receptor_col='gene_symbol', logfc_col='logFC', padj_col='padj',
                 deg_gene_col='gene', deg_logfc_col='logfoldchanges', deg_padj_col='pvals_adj'):
    """Run mismatch analysis for one comparison."""

    print(f"\n{'=' * 80}")
    print(f"Mismatch analysis: {comparison_name}")
    print(f"{'=' * 80}")

    # Upregulated receptors
    receptors_up = receptor_df[
        (receptor_df[logfc_col] > 0) &
        (receptor_df[padj_col] < 0.05)
    ].copy()
    print(f"Upregulated receptors (padj<0.05): {len(receptors_up)}")

    # Full gene lookup
    gene_lfc = deg_full_df.set_index(deg_gene_col)[deg_logfc_col].to_dict()
    gene_padj = deg_full_df.set_index(deg_gene_col)[deg_padj_col].to_dict()

    pair_results = []

    for _, row in receptors_up.iterrows():
        receptor = row[receptor_col]
        r_logfc = row[logfc_col]

        cognates = lr_exp[lr_exp['receptor'] == receptor][['ligand', 'n_db']]
        hc_cognates = cognates[cognates['n_db'] >= DB_THRESHOLD]

        if len(hc_cognates) == 0:
            continue

        lig_data = []
        for _, lr_row in hc_cognates.iterrows():
            lig = lr_row['ligand']
            ndb = lr_row['n_db']
            if lig in gene_lfc and np.isfinite(gene_lfc[lig]):
                lig_data.append({
                    'ligand': lig, 'logfc': gene_lfc[lig],
                    'padj': gene_padj.get(lig, 1.0), 'n_db': ndb
                })

        if len(lig_data) == 0:
            continue

        lig_df = pd.DataFrame(lig_data)
        n_hc = len(lig_df)
        sig_down = lig_df[(lig_df['logfc'] < 0) & (lig_df['padj'] < 0.05)]
        n_sig_down = len(sig_down)
        starvation = n_sig_down / n_hc

        for _, lr_row in sig_down.iterrows():
            r_c = np.clip(r_logfc, 0, 10) / 10.0
            l_c = np.clip(abs(lr_row['logfc']), 0, 10) / 10.0
            db_w = lr_row['n_db'] / 5.0
            composite = r_c * l_c * starvation * db_w

            pair_results.append({
                'receptor': receptor,
                'ligand': lr_row['ligand'],
                'r_logfc': round(r_logfc, 4),
                'l_logfc': round(lr_row['logfc'], 4),
                'l_padj': lr_row['padj'],
                'n_db': lr_row['n_db'],
                'db_weight': round(db_w, 2),
                'starvation': round(starvation, 3),
                'n_sig_down_hc': n_sig_down,
                'n_hc_ligands': n_hc,
                'composite': round(composite, 4),
            })

    pairs = pd.DataFrame(pair_results)
    if len(pairs) == 0:
        print("No mismatch pairs found!")
        return pairs

    pairs = pairs.sort_values('composite', ascending=False).reset_index(drop=True)
    pairs.index += 1
    pairs.index.name = 'rank'
    print(f"Total mismatch pairs: {len(pairs)}")

    # Display top 30
    cols = ['receptor', 'ligand', 'r_logfc', 'l_logfc', 'n_db',
            'starvation', 'n_sig_down_hc', 'n_hc_ligands', 'composite']
    print(f"\nTop 30:")
    print(pairs[cols].head(30).to_string())

    # Key pairs
    print(f"\nKey pairs:")
    for r, l in [('FGFR2','FGF10'),('FGFR2','FGF7'),('FGFR2','FGF1'),
                 ('BMPR2','BMP4'),('BMPR1A','BMP4'),('ACVR1','BMP6'),
                 ('FZD2','WNT9A'),('LRP6','WNT5B')]:
        m = pairs[(pairs['receptor']==r) & (pairs['ligand']==l)]
        if len(m) > 0:
            i = m.index[0]
            row = m.iloc[0]
            print(f"  {r}/{l:<10} rank={i:>4}/{len(pairs)}  comp={row['composite']:.4f}  "
                  f"starv={row['starvation']:.3f}  n_db={row['n_db']}")
        else:
            # Check why missing
            in_db = len(lr_exp[(lr_exp['receptor']==r) & (lr_exp['ligand']==l) & (lr_exp['n_db']>=DB_THRESHOLD)]) > 0
            r_up = receptor in receptors_up[receptor_col].values
            l_in = l in gene_lfc
            l_down = gene_lfc.get(l, 0) < 0 if l_in else False
            l_sig = gene_padj.get(l, 1) < 0.05 if l_in else False
            print(f"  {r}/{l:<10} NOT FOUND (in_db={in_db}, r_up={r_up}, "
                  f"l_in_deg={l_in}, l_down={l_down}, l_sig={l_sig})")

    return pairs


# --- A) Quiescent vs Activated ---
pairs_state = run_mismatch(
    rec_state, deg_wilcox, "Quiescent vs Activated (cell state)",
    receptor_col='gene_symbol', logfc_col='logFC', padj_col='padj',
    deg_gene_col='gene', deg_logfc_col='logfoldchanges', deg_padj_col='pvals_adj'
)

# --- B) Normal vs MI ---
pairs_cond = run_mismatch(
    rec_cond, deg_wilcox, "Normal vs MI (condition)",
    receptor_col='gene_symbol', logfc_col='logFC', padj_col='padj',
    deg_gene_col='gene', deg_logfc_col='logfoldchanges', deg_padj_col='pvals_adj'
)

# ============================================================================
# 4. Save human results
# ============================================================================
print("\n" + "=" * 80)
print("Saving results")
print("=" * 80)

if len(pairs_state) > 0:
    pairs_state.to_csv(OUTPUT_DIR / "human_lr_mismatch_state.csv")
    print(f"Saved: human_lr_mismatch_state.csv ({len(pairs_state)} pairs)")

if len(pairs_cond) > 0:
    pairs_cond.to_csv(OUTPUT_DIR / "human_lr_mismatch_condition.csv")
    print(f"Saved: human_lr_mismatch_condition.csv ({len(pairs_cond)} pairs)")
