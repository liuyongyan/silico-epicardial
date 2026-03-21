#!/usr/bin/env python3
"""
Ligand-Receptor Mismatch Analysis for Mouse Epicardial Data (Quaife-Ryan 2021)

Identifies "primed but starved" pairs: receptor↑ + ligand↓ in activated vs quiescent.

L-R Database: OmniPath merged from CellPhoneDB + CellTalkDB + Ramilowski2015 +
              Fantom5_LRdb + connectomeDB2020 (matches instruction 3 spec).
"""

import pandas as pd
import numpy as np
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent.parent.parent
RESULTS_DIR = PROJECT_DIR / "results" / "mouse"
OUTPUT_DIR = PROJECT_DIR / "results" / "mismatch"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================================
# 1. Load L-R database (5 sources merged, matching instruction 3)
# ============================================================================
print("=" * 70)
print("Step 1: Loading L-R database (CellPhoneDB+CellTalkDB+Ramilowski+Fantom5+connectomeDB)")
print("=" * 70)

import omnipath

lr = omnipath.interactions.PostTranslational.get(genesymbols=True, organism='human')

# Filter for the 5 databases specified in instruction 3
sources_filter = 'CellPhoneDB|CellTalkDB|Ramilowski2015|Fantom5_LRdb|connectomeDB2020'
lr_curated = lr[lr['sources'].str.contains(sources_filter, na=False)].copy()

# Extract raw pairs
raw_pairs = lr_curated[['source_genesymbol', 'target_genesymbol']].drop_duplicates()
raw_pairs.columns = ['ligand', 'receptor']
print(f"Raw curated pairs: {len(raw_pairs)}")

# ---- Handle complex receptor names (e.g., FZD2_LRP6 -> FZD2, LRP6) ----
expanded_rows = []
for _, row in raw_pairs.iterrows():
    lig = row['ligand']
    rec = row['receptor']

    # Split ligand complexes
    ligands = lig.split('_') if '_' in lig else [lig]
    # Split receptor complexes
    receptors = rec.split('_') if '_' in rec else [rec]

    for l in ligands:
        for r in receptors:
            expanded_rows.append({'ligand': l, 'receptor': r})

lr_expanded = pd.DataFrame(expanded_rows).drop_duplicates()
print(f"After expanding complexes: {len(lr_expanded)}")

# Convert human -> mouse gene symbols
def human_to_mouse(gene):
    if gene is None or pd.isna(gene):
        return None
    return gene[0].upper() + gene[1:].lower()

lr_pairs = lr_expanded.copy()
lr_pairs['ligand'] = lr_pairs['ligand'].apply(human_to_mouse)
lr_pairs['receptor'] = lr_pairs['receptor'].apply(human_to_mouse)
lr_pairs = lr_pairs.drop_duplicates()

print(f"Mouse L-R pairs: {len(lr_pairs)}")
print(f"Unique ligands: {lr_pairs['ligand'].nunique()}")
print(f"Unique receptors: {lr_pairs['receptor'].nunique()}")

# Verify key pairs
key_check = [('Fgf10', 'Fgfr2'), ('Bmp4', 'Bmpr2'), ('Bmp4', 'Bmpr1a'),
             ('Bmp6', 'Acvr1'), ('Wnt9a', 'Fzd2'), ('Wnt5b', 'Lrp6'),
             ('Fgf7', 'Fgfr2'), ('Fgf1', 'Fgfr2')]
print("\nKey pair verification:")
for lig, rec in key_check:
    found = len(lr_pairs[(lr_pairs['ligand'] == lig) & (lr_pairs['receptor'] == rec)]) > 0
    print(f"  {lig}/{rec}: {'FOUND' if found else 'MISSING'}")

# ============================================================================
# 2. Load mouse DEG data
# ============================================================================
print("\n" + "=" * 70)
print("Step 2: Loading mouse DEG data (Activated vs Quiescent)")
print("=" * 70)

deg_all = pd.read_csv(RESULTS_DIR / "quaife_ryan_de_activated_vs_quiescent.csv")
deg_all = deg_all.dropna(subset=['logfoldchanges'])
deg_all = deg_all[np.isfinite(deg_all['logfoldchanges'])]
print(f"Genes with valid logFC: {len(deg_all)}")

receptor_ranks = pd.read_csv(RESULTS_DIR / "receptor_rankings_by_logfc.csv")
receptors_up = receptor_ranks[
    (receptor_ranks['logfoldchanges'] > 0) &
    (receptor_ranks['pvals_adj'] < 0.05)
].copy()
print(f"Significantly upregulated receptors: {len(receptors_up)}")

# ============================================================================
# 3. Find mismatch pairs
# ============================================================================
print("\n" + "=" * 70)
print("Step 3: Finding 'Primed But Starved' mismatch pairs")
print("=" * 70)

gene_logfc = deg_all.set_index('names')['logfoldchanges'].to_dict()
gene_padj = deg_all.set_index('names')['pvals_adj'].to_dict()

mismatched = []

for _, row in receptors_up.iterrows():
    receptor = row['names']
    receptor_logfc = row['logfoldchanges']
    receptor_padj = row['pvals_adj']
    receptor_rank = row['rank']

    cognate_ligands = lr_pairs[lr_pairs['receptor'] == receptor]['ligand'].unique()

    for ligand in cognate_ligands:
        if ligand in gene_logfc:
            ligand_logfc = gene_logfc[ligand]
            ligand_padj = gene_padj.get(ligand, 1.0)

            if not np.isfinite(ligand_logfc):
                continue

            if ligand_logfc < 0:
                mismatch_score = receptor_logfc - ligand_logfc
                mismatched.append({
                    'receptor': receptor,
                    'ligand': ligand,
                    'receptor_logfc': round(receptor_logfc, 2),
                    'receptor_padj': receptor_padj,
                    'receptor_rank_among_receptors': int(receptor_rank),
                    'ligand_logfc': round(ligand_logfc, 2),
                    'ligand_padj': ligand_padj,
                    'mismatch_score': round(mismatch_score, 2),
                })

mismatch_df = pd.DataFrame(mismatched)
mismatch_df = mismatch_df.sort_values('mismatch_score', ascending=False).reset_index(drop=True)
mismatch_df.index = mismatch_df.index + 1
mismatch_df.index.name = 'rank'

print(f"Total mismatch pairs: {len(mismatch_df)}")
print(f"Unique receptors: {mismatch_df['receptor'].nunique()}")
print(f"Unique ligands: {mismatch_df['ligand'].nunique()}")

# ============================================================================
# 4. Annotate pathways
# ============================================================================
pathway_map = {
    'Fgfr1': 'FGF', 'Fgfr2': 'FGF', 'Fgfr3': 'FGF', 'Fgfrl1': 'FGF',
    'Bmpr1a': 'BMP', 'Bmpr1b': 'BMP', 'Bmpr2': 'BMP',
    'Acvr1': 'BMP/Activin', 'Acvr1b': 'Activin', 'Acvr2a': 'Activin', 'Acvr2b': 'Activin',
    'Tgfbr1': 'TGF-β', 'Tgfbr2': 'TGF-β', 'Tgfbr3': 'TGF-β',
    'Fzd1': 'Wnt', 'Fzd2': 'Wnt', 'Fzd3': 'Wnt', 'Fzd4': 'Wnt',
    'Fzd5': 'Wnt', 'Fzd6': 'Wnt', 'Fzd7': 'Wnt', 'Fzd8': 'Wnt',
    'Lrp5': 'Wnt', 'Lrp6': 'Wnt', 'Ryk': 'Wnt', 'Ror2': 'Wnt',
    'Pdgfra': 'PDGF', 'Pdgfrb': 'PDGF', 'Pdgfrl': 'PDGF',
    'Egfr': 'EGF', 'Erbb2': 'EGF', 'Erbb3': 'EGF', 'Erbb4': 'EGF',
    'Kdr': 'VEGF', 'Flt1': 'VEGF', 'Flt4': 'VEGF',
    'Notch1': 'Notch', 'Notch2': 'Notch', 'Notch3': 'Notch', 'Notch4': 'Notch',
    'Met': 'HGF', 'Igf1r': 'IGF', 'Igf2r': 'IGF',
    'Nrp1': 'Semaphorin/VEGF', 'Nrp2': 'Semaphorin/VEGF',
    'Plxna1': 'Semaphorin', 'Plxna2': 'Semaphorin', 'Plxna4': 'Semaphorin',
    'Plxnb1': 'Semaphorin', 'Plxnb2': 'Semaphorin',
}
mismatch_df['pathway'] = mismatch_df['receptor'].map(pathway_map).fillna('Other')

# ============================================================================
# 5. Display top results
# ============================================================================
print("\n" + "=" * 70)
print("Step 4: Top 40 'Primed But Starved' pairs")
print("=" * 70)

display_cols = ['receptor', 'ligand', 'receptor_logfc', 'ligand_logfc',
                'mismatch_score', 'pathway', 'receptor_rank_among_receptors']
print(mismatch_df[display_cols].head(40).to_string())

# ============================================================================
# 6. Key pairs from instruction 3
# ============================================================================
print("\n" + "=" * 70)
print("Step 5: Instruction 3 key pairs — actual vs claimed ranks")
print("=" * 70)

key_pairs = [
    ('Bmpr2', 'Bmp4', 1),
    ('Bmpr1a', 'Bmp4', 2),
    ('Acvr1', 'Bmp6', 3),
    ('Fzd2', 'Wnt9a', 4),
    ('Lrp6', 'Wnt5b', 5),
    ('Fgfr2', 'Fgf10', 6),
    ('Fgfr2', 'Fgf7', 7),
    ('Fgfr2', 'Fgf1', 8),
]

print(f"\n{'Pair':<20} {'Inst3':>6} {'Actual':>8} {'Overall':>10} {'R_logFC':>9} {'L_logFC':>9} {'Score':>8}")
print("-" * 76)

for receptor, ligand, inst3_rank in key_pairs:
    match = mismatch_df[(mismatch_df['receptor'] == receptor) &
                         (mismatch_df['ligand'] == ligand)]
    if len(match) > 0:
        idx = match.index[0]
        row = match.iloc[0]
        total = len(mismatch_df)
        print(f"  {receptor}/{ligand:<13} {inst3_rank:>5} {idx:>7}/{total} "
              f"{row['receptor_logfc']:>9.2f} {row['ligand_logfc']:>9.2f} {row['mismatch_score']:>8.2f}")
    else:
        in_db = len(lr_pairs[(lr_pairs['receptor'] == receptor) & (lr_pairs['ligand'] == ligand)]) > 0
        print(f"  {receptor}/{ligand:<13} {inst3_rank:>5} {'N/A':>8}   (in_LR_db={in_db})")

# ============================================================================
# 7. Best pair per pathway
# ============================================================================
print("\n" + "=" * 70)
print("Step 6: Best mismatch pair per signaling pathway")
print("=" * 70)

for pathway in ['BMP', 'BMP/Activin', 'FGF', 'Wnt', 'TGF-β', 'EGF', 'PDGF',
                'Notch', 'HGF', 'IGF', 'Semaphorin']:
    pw = mismatch_df[mismatch_df['pathway'] == pathway]
    if len(pw) > 0:
        best = pw.iloc[0]
        rank = pw.index[0]
        print(f"  {pathway:<15} #{rank:<5} {best['receptor']}/{best['ligand']:<20} "
              f"R={best['receptor_logfc']:>7.2f}  L={best['ligand_logfc']:>7.2f}  "
              f"score={best['mismatch_score']:>7.2f}  ({len(pw)} pairs)")

# ============================================================================
# 8. Save
# ============================================================================
print("\n" + "=" * 70)
print("Step 7: Saving results")
print("=" * 70)

mismatch_df.to_csv(OUTPUT_DIR / "mouse_lr_mismatch_all.csv")
mismatch_df.head(50).to_csv(OUTPUT_DIR / "mouse_lr_mismatch_top50.csv")
lr_pairs.to_csv(OUTPUT_DIR / "curated_lr_pairs_mouse.csv", index=False)

print(f"Saved: mouse_lr_mismatch_all.csv ({len(mismatch_df)} pairs)")
print(f"Saved: mouse_lr_mismatch_top50.csv")
print(f"Saved: curated_lr_pairs_mouse.csv ({len(lr_pairs)} pairs)")
print(f"\nL-R database: {len(lr_pairs)} pairs from {sources_filter}")
