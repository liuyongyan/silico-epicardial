#!/usr/bin/env python3
"""
Cross-species L-R mismatch comparison: Mouse vs Human

Human data: PERIHEART epicardial cells (Wilcoxon DEG, quiescent vs activated)
Mouse data: Quaife-Ryan epicardial cells (Wilcoxon DEG, activated vs quiescent)

Uses two significance tiers:
  - Strict: both receptor padj<0.05 AND ligand padj<0.05
  - Relaxed: receptor padj<0.05, ligand logFC<0 (any padj)
"""

import pandas as pd
import numpy as np
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent.parent.parent
OUTPUT_DIR = PROJECT_DIR / "results" / "mismatch"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================================
# 1. Load L-R database
# ============================================================================
print("=" * 80)
print("Step 1: Loading L-R database")
print("=" * 80)

import omnipath

lr_raw = omnipath.interactions.PostTranslational.get(genesymbols=True, organism='human')
sources_filter = 'CellPhoneDB|CellTalkDB|Ramilowski2015|Fantom5_LRdb|connectomeDB2020'
lr_curated = lr_raw[lr_raw['sources'].str.contains(sources_filter, na=False)].copy()

target_dbs = ['CellPhoneDB', 'CellTalkDB', 'Ramilowski2015', 'Fantom5_LRdb', 'connectomeDB2020']

def count_dbs(s):
    return sum(1 for db in target_dbs if db in s)

lr_curated['n_db'] = lr_curated['sources'].apply(count_dbs)

raw_pairs = lr_curated[['source_genesymbol', 'target_genesymbol', 'n_db']].copy()
raw_pairs.columns = ['ligand', 'receptor', 'n_db']

expanded = []
for _, row in raw_pairs.iterrows():
    for l in (row['ligand'].split('_') if '_' in row['ligand'] else [row['ligand']]):
        for r in (row['receptor'].split('_') if '_' in row['receptor'] else [row['receptor']]):
            expanded.append({'ligand': l, 'receptor': r, 'n_db': row['n_db']})

lr_db = pd.DataFrame(expanded).groupby(['ligand', 'receptor'])['n_db'].max().reset_index()
lr_hc = lr_db[lr_db['n_db'] >= 3]  # high-confidence only
print(f"L-R pairs (n_db>=3): {len(lr_hc)}")

# ============================================================================
# 2. Load DEG data for both species
# ============================================================================
print("\n" + "=" * 80)
print("Step 2: Loading DEG data")
print("=" * 80)

# --- Human (Wilcoxon, quiescent vs activated) ---
h_deg = pd.read_csv(PROJECT_DIR / "results" / "deg" / "deg_results_wilcoxon.csv")
h_deg = h_deg.dropna(subset=['logfoldchanges'])
h_deg = h_deg[np.isfinite(h_deg['logfoldchanges'])]
h_lfc = h_deg.set_index('gene')['logfoldchanges'].to_dict()
h_padj = h_deg.set_index('gene')['pvals_adj'].to_dict()
print(f"Human genes: {len(h_deg)}")

# --- Mouse (Wilcoxon, activated vs quiescent) ---
m_deg = pd.read_csv(PROJECT_DIR / "results" / "mouse" / "quaife_ryan_de_activated_vs_quiescent.csv")
m_deg = m_deg.dropna(subset=['logfoldchanges'])
m_deg_finite = m_deg[np.isfinite(m_deg['logfoldchanges'])]
m_lfc = m_deg_finite.set_index('names')['logfoldchanges'].to_dict()
m_padj = m_deg_finite.set_index('names')['pvals_adj'].to_dict()
print(f"Mouse genes: {len(m_deg_finite)}")

# Mouse gene -> human ortholog (simple case conversion for most genes)
def mouse_to_human(g):
    return g.upper() if g else None

def human_to_mouse(g):
    return g[0].upper() + g[1:].lower() if g else None

# ============================================================================
# 3. For each high-confidence L-R pair, check pattern in both species
# ============================================================================
print("\n" + "=" * 80)
print("Step 3: Cross-species mismatch comparison")
print("=" * 80)

results = []

for _, lr_row in lr_hc.iterrows():
    h_receptor = lr_row['receptor']  # human gene symbol
    h_ligand = lr_row['ligand']
    n_db = lr_row['n_db']

    m_receptor = human_to_mouse(h_receptor)
    m_ligand = human_to_mouse(h_ligand)

    # Human data
    h_r_lfc = h_lfc.get(h_receptor, None)
    h_r_padj = h_padj.get(h_receptor, None)
    h_l_lfc = h_lfc.get(h_ligand, None)
    h_l_padj = h_padj.get(h_ligand, None)

    # Mouse data
    m_r_lfc = m_lfc.get(m_receptor, None)
    m_r_padj = m_padj.get(m_receptor, None)
    m_l_lfc = m_lfc.get(m_ligand, None)
    m_l_padj = m_padj.get(m_ligand, None)

    # Check "primed but starved" pattern
    # Receptor UP + Ligand DOWN
    def check_pattern(r_lfc, r_padj, l_lfc, l_padj, strict=True):
        if r_lfc is None or l_lfc is None:
            return 'no_data'
        if not np.isfinite(r_lfc) or not np.isfinite(l_lfc):
            return 'no_data'
        r_up = r_lfc > 0 and (r_padj < 0.05 if strict else True)
        l_down = l_lfc < 0 and (l_padj < 0.05 if strict else True)
        if r_up and l_down:
            return 'yes'
        elif r_lfc > 0 and l_lfc < 0:
            return 'trend'  # direction correct but not significant
        else:
            return 'no'

    m_strict = check_pattern(m_r_lfc, m_r_padj, m_l_lfc, m_l_padj, strict=True)
    m_relaxed = check_pattern(m_r_lfc, m_r_padj, m_l_lfc, m_l_padj, strict=False)
    h_strict = check_pattern(h_r_lfc, h_r_padj, h_l_lfc, h_l_padj, strict=True)
    h_relaxed = check_pattern(h_r_lfc, h_r_padj, h_l_lfc, h_l_padj, strict=False)

    results.append({
        'receptor': h_receptor,
        'ligand': h_ligand,
        'n_db': n_db,
        # Mouse
        'm_r_lfc': round(m_r_lfc, 4) if m_r_lfc is not None and np.isfinite(m_r_lfc) else None,
        'm_r_padj': m_r_padj,
        'm_l_lfc': round(m_l_lfc, 4) if m_l_lfc is not None and np.isfinite(m_l_lfc) else None,
        'm_l_padj': m_l_padj,
        'm_strict': m_strict,
        'm_relaxed': m_relaxed,
        # Human
        'h_r_lfc': round(h_r_lfc, 4) if h_r_lfc is not None else None,
        'h_r_padj': h_r_padj,
        'h_l_lfc': round(h_l_lfc, 4) if h_l_lfc is not None else None,
        'h_l_padj': h_l_padj,
        'h_strict': h_strict,
        'h_relaxed': h_relaxed,
    })

df = pd.DataFrame(results)

# ============================================================================
# 4. Categorize conservation
# ============================================================================
def conservation_category(row):
    if row['m_strict'] == 'yes' and row['h_strict'] == 'yes':
        return 'CONSERVED_strict'
    elif row['m_relaxed'] == 'yes' and row['h_relaxed'] == 'yes':
        return 'CONSERVED_relaxed'
    elif (row['m_strict'] == 'yes' or row['m_relaxed'] == 'yes') and row['h_relaxed'] in ('yes', 'trend'):
        return 'CONSERVED_trend'
    elif row['m_strict'] == 'yes' or row['m_relaxed'] == 'yes':
        return 'mouse_only'
    elif row['h_strict'] == 'yes' or row['h_relaxed'] == 'yes':
        return 'human_only'
    else:
        return 'neither'

df['conservation'] = df.apply(conservation_category, axis=1)

# ============================================================================
# 5. Display results
# ============================================================================

# Summary
print("\nConservation summary:")
print(df['conservation'].value_counts().to_string())

# --- CONSERVED pairs (both species show pattern) ---
for cat in ['CONSERVED_strict', 'CONSERVED_relaxed', 'CONSERVED_trend']:
    conserved = df[df['conservation'] == cat].sort_values('n_db', ascending=False)
    if len(conserved) > 0:
        print(f"\n{'=' * 80}")
        print(f"{cat}: {len(conserved)} pairs")
        print(f"{'=' * 80}")
        cols = ['receptor', 'ligand', 'n_db',
                'm_r_lfc', 'm_l_lfc', 'm_strict',
                'h_r_lfc', 'h_l_lfc', 'h_strict']
        print(conserved[cols].head(40).to_string(index=False))

# --- Key pairs ---
print(f"\n{'=' * 80}")
print("Key pairs from instruction 3")
print(f"{'=' * 80}")
print(f"\n{'Pair':<18} {'n_db':>4}  {'Mouse_R':>8} {'Mouse_L':>8} {'M_pat':>8}  "
      f"{'Human_R':>8} {'Human_L':>8} {'H_pat':>8}  {'Conserv':<20}")
print("-" * 110)

for r, l in [('FGFR2','FGF10'),('FGFR2','FGF7'),('FGFR2','FGF1'),
             ('BMPR2','BMP4'),('BMPR1A','BMP4'),('ACVR1','BMP6'),
             ('FZD2','WNT9A'),('LRP6','WNT5B'),
             ('FGFR2','FGF2'),('BMPR2','BMP6'),('BMPR2','BMP2')]:
    m = df[(df['receptor']==r) & (df['ligand']==l)]
    if len(m) > 0:
        row = m.iloc[0]
        mr = f"{row['m_r_lfc']:>8.2f}" if row['m_r_lfc'] is not None else '     N/A'
        ml = f"{row['m_l_lfc']:>8.2f}" if row['m_l_lfc'] is not None else '     N/A'
        hr = f"{row['h_r_lfc']:>8.2f}" if row['h_r_lfc'] is not None else '     N/A'
        hl = f"{row['h_l_lfc']:>8.2f}" if row['h_l_lfc'] is not None else '     N/A'
        print(f"  {r}/{l:<11} {row['n_db']:>4}  {mr} {ml} {row['m_relaxed']:>8}  "
              f"{hr} {hl} {row['h_relaxed']:>8}  {row['conservation']:<20}")
    else:
        print(f"  {r}/{l:<11}   not in HC database")

# ============================================================================
# 6. Mouse-only pairs (potential novel targets if human data were better)
# ============================================================================
mouse_only = df[df['conservation'] == 'mouse_only'].sort_values('n_db', ascending=False)
print(f"\n{'=' * 80}")
print(f"Mouse-only 'primed but starved' (top 20 by n_db): {len(mouse_only)} pairs")
print(f"{'=' * 80}")
cols = ['receptor', 'ligand', 'n_db', 'm_r_lfc', 'm_l_lfc', 'm_strict',
        'h_r_lfc', 'h_l_lfc', 'h_relaxed']
print(mouse_only[cols].head(20).to_string(index=False))

# ============================================================================
# 7. Save
# ============================================================================
print(f"\n{'=' * 80}")
print("Saving")
print(f"{'=' * 80}")

df.to_csv(OUTPUT_DIR / "cross_species_lr_mismatch.csv", index=False)
print(f"Saved: cross_species_lr_mismatch.csv ({len(df)} pairs)")

# Save conserved pairs
conserved_all = df[df['conservation'].str.startswith('CONSERVED')]
conserved_all.to_csv(OUTPUT_DIR / "cross_species_conserved.csv", index=False)
print(f"Saved: cross_species_conserved.csv ({len(conserved_all)} pairs)")
