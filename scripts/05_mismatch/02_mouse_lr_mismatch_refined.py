#!/usr/bin/env python3
"""
Refined L-R Mismatch Analysis with multi-dimensional scoring.

Scoring dimensions:
  1. Receptor logFC (capped ±10, normalized)
  2. Ligand logFC (capped ±10, abs, normalized) — must be sig down (padj<0.05)
  3. Starvation ratio: fraction of HIGH-CONFIDENCE cognate ligands that are sig down
  4. Database consensus: n_db/5 (how many of 5 curated DBs support the pair)

composite = r_norm × l_norm × starvation_ratio × db_weight
"""

import pandas as pd
import numpy as np
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent.parent.parent
RESULTS_DIR = PROJECT_DIR / "results" / "mouse"
OUTPUT_DIR = PROJECT_DIR / "results" / "mismatch"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================================
# 1. Build L-R database with per-pair database consensus count
# ============================================================================
print("=" * 80)
print("Step 1: Building L-R database with database consensus scores")
print("=" * 80)

import omnipath

lr_raw = omnipath.interactions.PostTranslational.get(genesymbols=True, organism='human')
sources_filter = 'CellPhoneDB|CellTalkDB|Ramilowski2015|Fantom5_LRdb|connectomeDB2020'
lr_curated = lr_raw[lr_raw['sources'].str.contains(sources_filter, na=False)].copy()

target_dbs = ['CellPhoneDB', 'CellTalkDB', 'Ramilowski2015', 'Fantom5_LRdb', 'connectomeDB2020']

def count_dbs(sources_str):
    return sum(1 for db in target_dbs if db in sources_str)

lr_curated['n_db'] = lr_curated['sources'].apply(count_dbs)

# Extract ligand->receptor pairs with n_db
raw_pairs = lr_curated[['source_genesymbol', 'target_genesymbol', 'n_db']].copy()
raw_pairs.columns = ['ligand', 'receptor', 'n_db']

# Expand complexes
expanded = []
for _, row in raw_pairs.iterrows():
    ligs = row['ligand'].split('_') if '_' in row['ligand'] else [row['ligand']]
    recs = row['receptor'].split('_') if '_' in row['receptor'] else [row['receptor']]
    for l in ligs:
        for r in recs:
            expanded.append({'ligand': l, 'receptor': r, 'n_db': row['n_db']})

lr_exp = pd.DataFrame(expanded)
# Keep max n_db per pair
lr_exp = lr_exp.groupby(['ligand', 'receptor'])['n_db'].max().reset_index()

# Convert to mouse
def h2m(g):
    return g[0].upper() + g[1:].lower() if g else None

lr_exp['ligand'] = lr_exp['ligand'].apply(h2m)
lr_exp['receptor'] = lr_exp['receptor'].apply(h2m)
lr_exp = lr_exp.drop_duplicates()

print(f"L-R pairs with consensus scores: {len(lr_exp)}")
print(f"Distribution of n_db: {lr_exp['n_db'].value_counts().sort_index().to_dict()}")

# ============================================================================
# 2. Load DEG data
# ============================================================================
print("\n" + "=" * 80)
print("Step 2: Loading DEG data")
print("=" * 80)

deg = pd.read_csv(RESULTS_DIR / "quaife_ryan_de_activated_vs_quiescent.csv")
deg = deg.dropna(subset=['logfoldchanges'])
gene_lfc = deg.set_index('names')['logfoldchanges'].to_dict()
gene_padj = deg.set_index('names')['pvals_adj'].to_dict()

rec = pd.read_csv(RESULTS_DIR / "receptor_rankings_by_logfc.csv")
receptors_up = rec[(rec['logfoldchanges'] > 0) & (rec['pvals_adj'] < 0.05)]
print(f"Upregulated receptors: {len(receptors_up)}")

# ============================================================================
# 3. Compute starvation ratio using ONLY high-confidence ligands (n_db >= 3)
# ============================================================================
print("\n" + "=" * 80)
print("Step 3: Computing starvation ratio (high-confidence ligands, n_db >= 3)")
print("=" * 80)

DB_THRESHOLD = 3  # only count ligands supported by >= 3 databases

pair_results = []

for _, row in receptors_up.iterrows():
    receptor = row['names']
    r_logfc = row['logfoldchanges']
    r_padj = row['pvals_adj']

    # Get cognate ligands with their database support
    cognates = lr_exp[lr_exp['receptor'] == receptor][['ligand', 'n_db']]

    # Filter for high-confidence pairs
    hc_cognates = cognates[cognates['n_db'] >= DB_THRESHOLD]

    if len(hc_cognates) == 0:
        continue

    # Check expression of high-confidence ligands
    lig_data = []
    for _, lr_row in hc_cognates.iterrows():
        lig = lr_row['ligand']
        ndb = lr_row['n_db']
        if lig in gene_lfc and np.isfinite(gene_lfc[lig]):
            lig_data.append({
                'ligand': lig,
                'logfc': gene_lfc[lig],
                'padj': gene_padj.get(lig, 1.0),
                'n_db': ndb
            })

    if len(lig_data) == 0:
        continue

    lig_df = pd.DataFrame(lig_data)
    n_hc = len(lig_df)
    sig_down = lig_df[(lig_df['logfc'] < 0) & (lig_df['padj'] < 0.05)]
    n_sig_down = len(sig_down)
    starvation = n_sig_down / n_hc

    for _, lr_row in sig_down.iterrows():
        r_c = np.clip(r_logfc, 0, 10)
        l_c = np.clip(abs(lr_row['logfc']), 0, 10)
        db_weight = lr_row['n_db'] / 5.0

        composite = (r_c / 10) * (l_c / 10) * starvation * db_weight

        pair_results.append({
            'receptor': receptor,
            'ligand': lr_row['ligand'],
            'r_logfc': round(r_logfc, 2),
            'l_logfc': round(lr_row['logfc'], 2),
            'l_padj': lr_row['padj'],
            'n_db': lr_row['n_db'],
            'db_weight': round(db_weight, 2),
            'starvation': round(starvation, 3),
            'n_sig_down_hc': n_sig_down,
            'n_hc_ligands': n_hc,
            'composite': round(composite, 4),
        })

pairs = pd.DataFrame(pair_results)
pairs = pairs.sort_values('composite', ascending=False).reset_index(drop=True)
pairs.index += 1
pairs.index.name = 'rank'

print(f"Total mismatch pairs: {len(pairs)}")

# ============================================================================
# 4. Pathway annotation
# ============================================================================
pathway_map = {
    'Fgfr1': 'FGF', 'Fgfr2': 'FGF', 'Fgfr3': 'FGF', 'Fgfrl1': 'FGF',
    'Bmpr1a': 'BMP', 'Bmpr1b': 'BMP', 'Bmpr2': 'BMP',
    'Acvr1': 'BMP/Activin', 'Acvr1b': 'Activin', 'Acvr2a': 'Activin',
    'Tgfbr1': 'TGF-β', 'Tgfbr2': 'TGF-β', 'Tgfbr3': 'TGF-β',
    'Fzd1': 'Wnt', 'Fzd2': 'Wnt', 'Fzd3': 'Wnt', 'Fzd4': 'Wnt',
    'Fzd5': 'Wnt', 'Fzd6': 'Wnt', 'Fzd7': 'Wnt', 'Fzd8': 'Wnt',
    'Lrp5': 'Wnt', 'Lrp6': 'Wnt', 'Ryk': 'Wnt', 'Ror2': 'Wnt',
    'Pdgfra': 'PDGF', 'Pdgfrb': 'PDGF', 'Pdgfrl': 'PDGF',
    'Egfr': 'EGF', 'Erbb2': 'EGF', 'Erbb3': 'EGF', 'Erbb4': 'EGF',
    'Kdr': 'VEGF', 'Flt1': 'VEGF', 'Flt4': 'VEGF',
    'Notch1': 'Notch', 'Notch2': 'Notch', 'Notch3': 'Notch', 'Notch4': 'Notch',
    'Met': 'HGF', 'Igf1r': 'IGF', 'Igf2r': 'IGF',
    'Nrp1': 'Sema/VEGF', 'Nrp2': 'Sema/VEGF',
    'Plxna1': 'Sema', 'Plxna2': 'Sema', 'Plxna4': 'Sema',
}
pairs['pathway'] = pairs['receptor'].map(pathway_map).fillna('Other')

# ============================================================================
# 5. Display results
# ============================================================================
print("\n" + "=" * 80)
print("Top 40 pairs (composite = r_norm × l_norm × starvation × db_weight)")
print("=" * 80)

cols = ['receptor', 'ligand', 'r_logfc', 'l_logfc', 'n_db',
        'starvation', 'n_sig_down_hc', 'n_hc_ligands', 'composite', 'pathway']
print(pairs[cols].head(40).to_string())

# Key pairs
print("\n" + "=" * 80)
print("Key pairs from instruction 3")
print("=" * 80)
print(f"\n{'Pair':<20} {'Rank':>8} {'Composite':>10} {'Starv':>6} {'n_db':>5} {'HC_down/HC_tot':>15}")
print("-" * 70)
for r, l in [('Fgfr2','Fgf10'),('Fgfr2','Fgf7'),('Fgfr2','Fgf1'),
             ('Bmpr2','Bmp4'),('Bmpr1a','Bmp4'),('Acvr1','Bmp6'),
             ('Fzd2','Wnt9a'),('Lrp6','Wnt5b'),
             ('Itgb1','Spon2'),('Cd44','Lpl'),('Sdc1','Lpl')]:
    m = pairs[(pairs['receptor']==r) & (pairs['ligand']==l)]
    if len(m) > 0:
        i = m.index[0]
        row = m.iloc[0]
        print(f"  {r}/{l:<13} {i:>6}/{len(pairs)} {row['composite']:>10.4f} "
              f"{row['starvation']:>5.3f} {row['n_db']:>5} "
              f"{row['n_sig_down_hc']:>6}/{row['n_hc_ligands']}")
    else:
        print(f"  {r}/{l:<13} {'FILTERED OUT':>8}")

# Best per pathway
print("\n" + "=" * 80)
print("Best pair per pathway")
print("=" * 80)
for pw in ['BMP', 'BMP/Activin', 'FGF', 'Wnt', 'TGF-β', 'EGF', 'Notch', 'IGF', 'Sema']:
    sub = pairs[pairs['pathway'] == pw]
    if len(sub) > 0:
        best = sub.iloc[0]
        rank = sub.index[0]
        print(f"  {pw:<15} #{rank:<4} {best['receptor']}/{best['ligand']:<18} "
              f"comp={best['composite']:.4f}  starv={best['starvation']:.3f}  "
              f"n_db={best['n_db']}  ({len(sub)} pairs)")

# ============================================================================
# 6. Save
# ============================================================================
pairs.to_csv(OUTPUT_DIR / "mouse_lr_mismatch_refined.csv")
print(f"\nSaved: mouse_lr_mismatch_refined.csv ({len(pairs)} pairs)")
