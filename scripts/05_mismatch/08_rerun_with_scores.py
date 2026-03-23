#!/usr/bin/env python3
"""
Re-run full mismatch pipeline using Wilcoxon scores instead of logFC.

Wilcoxon score (z-normalized U statistic) is robust to near-zero expression
artifacts that inflate logFC. This replaces logFC in:
  - Step 2: Mouse mismatch composite
  - Step 3: Cross-species conservation ranking
  - Step 7: Final prioritization

All other steps (L-R database, starvation ratio, n_db filter, family filter,
DGIdb, UniProt, PubMed) remain unchanged.
"""

import pandas as pd
import numpy as np
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent.parent.parent
MISMATCH_DIR = PROJECT_DIR / "results" / "mismatch"
OUTPUT_DIR = MISMATCH_DIR

# ============================================================================
# 1. Load all data
# ============================================================================
print("=" * 80)
print("Step 1: Loading data")
print("=" * 80)

# Mouse DEG (has both scores and logfoldchanges)
m_deg = pd.read_csv(PROJECT_DIR / "results" / "mouse" / "quaife_ryan_de_activated_vs_quiescent.csv")
m_deg = m_deg.dropna(subset=['logfoldchanges', 'scores'])
print(f"Mouse DEG: {len(m_deg)} genes")

# Mouse receptor rankings
m_rec = pd.read_csv(PROJECT_DIR / "results" / "mouse" / "receptor_rankings_by_logfc.csv")
print(f"Mouse receptors: {len(m_rec)}")

# Human DEG
h_deg = pd.read_csv(PROJECT_DIR / "results" / "deg" / "deg_results_wilcoxon.csv")
h_deg = h_deg.dropna(subset=['logfoldchanges', 'scores'])
print(f"Human DEG: {len(h_deg)} genes")

# L-R database (already built, with n_db)
# Helper: human to mouse gene name
def h2m(g):
    return g[0].upper() + g[1:].lower() if g else None

# Use cross-species file as L-R database source (has n_db, human gene names)
cross_saved = pd.read_csv(MISMATCH_DIR / "cross_species_lr_mismatch.csv")
lr_db = cross_saved[['receptor', 'ligand', 'n_db']].drop_duplicates()
lr_db = lr_db.rename(columns={'receptor': 'receptor_human', 'ligand': 'ligand_human'})
lr_db['receptor'] = lr_db['receptor_human'].apply(h2m)
lr_db['ligand'] = lr_db['ligand_human'].apply(h2m)
lr_hc = lr_db[lr_db['n_db'] >= 3]
print(f"L-R pairs (n_db>=3): {len(lr_hc)}")

# Gene scores (DGIdb + UniProt + PubMed, already computed)
gene_scores = pd.read_csv(MISMATCH_DIR / "gene_scores_corrected.csv", index_col=0)
deliver_data = pd.read_csv(MISMATCH_DIR / "uniprot_deliverability.csv", index_col=0)
print(f"Gene scores: {len(gene_scores)}")

# ============================================================================
# 2. Mouse mismatch with Wilcoxon scores
# ============================================================================
print("\n" + "=" * 80)
print("Step 2: Mouse mismatch identification (using Wilcoxon scores)")
print("=" * 80)

# Upregulated receptors: score > 0 AND padj < 0.05
# (score > 0 = activated > quiescent in Wilcoxon rank-sum)
receptors_up = m_rec[(m_rec['scores'] > 0) & (m_rec['pvals_adj'] < 0.05)].copy()
print(f"Upregulated receptors (score>0, padj<0.05): {len(receptors_up)}")

# Build gene lookups using SCORES (not logFC)
m_score = m_deg.set_index('names')['scores'].to_dict()
m_padj = m_deg.set_index('names')['pvals_adj'].to_dict()
m_logfc = m_deg.set_index('names')['logfoldchanges'].to_dict()  # keep for reference

# Convert mouse gene to human for L-R lookup
def h2m(g):
    return g[0].upper() + g[1:].lower() if g else None

DB_THRESHOLD = 3

pair_results = []

for _, row in receptors_up.iterrows():
    receptor = row['names']
    r_score = row['scores']

    # Get HC cognate ligands (lr_hc has mouse gene names in 'ligand'/'receptor')
    cognates = lr_hc[lr_hc['receptor'] == receptor][['ligand', 'ligand_human', 'n_db']].copy()
    cognates = cognates.rename(columns={'ligand': 'ligand_mouse'})

    if len(cognates) == 0:
        continue

    # Check which HC ligands are significantly downregulated (score < 0, padj < 0.05)
    lig_data = []
    for _, lr_row in cognates.iterrows():
        lig_m = lr_row['ligand_mouse']
        ndb = lr_row['n_db']
        if lig_m in m_score:
            lig_data.append({
                'ligand': lig_m,
                'ligand_human': lr_row['ligand_human'],
                'score': m_score[lig_m],
                'padj': m_padj.get(lig_m, 1.0),
                'logfc': m_logfc.get(lig_m, 0),
                'n_db': ndb,
            })

    if len(lig_data) == 0:
        continue

    lig_df = pd.DataFrame(lig_data)
    n_hc = len(lig_df)
    sig_down = lig_df[(lig_df['score'] < 0) & (lig_df['padj'] < 0.05)]
    n_sig_down = len(sig_down)
    starvation = n_sig_down / n_hc

    for _, lr_row in sig_down.iterrows():
        # Normalize scores to [0, 1] range
        # Max receptor score ~190, max |ligand score| ~270
        r_norm = min(r_score / 200.0, 1.0)
        l_norm = min(abs(lr_row['score']) / 200.0, 1.0)
        db_w = lr_row['n_db'] / 5.0

        composite = r_norm * l_norm * starvation * db_w

        pair_results.append({
            'receptor': receptor,
            'ligand': lr_row['ligand'],
            'r_score': round(r_score, 2),
            'r_logfc': round(row['logfoldchanges'], 2),
            'l_score': round(lr_row['score'], 2),
            'l_logfc': round(lr_row['logfc'], 2),
            'l_padj': lr_row['padj'],
            'n_db': lr_row['n_db'],
            'starvation': round(starvation, 3),
            'n_sig_down_hc': n_sig_down,
            'n_hc_ligands': n_hc,
            'composite': round(composite, 4),
        })

mouse_pairs = pd.DataFrame(pair_results)
mouse_pairs = mouse_pairs.sort_values('composite', ascending=False).reset_index(drop=True)
mouse_pairs.index += 1
mouse_pairs.index.name = 'rank'

# Add pathway
pathway_map = {
    'Fgfr1':'FGF','Fgfr2':'FGF','Fgfr3':'FGF','Fgfrl1':'FGF',
    'Bmpr1a':'BMP','Bmpr1b':'BMP','Bmpr2':'BMP',
    'Acvr1':'BMP/Activin','Acvr1b':'Activin',
    'Tgfbr1':'TGF-β','Tgfbr2':'TGF-β',
    'Fzd1':'Wnt','Fzd2':'Wnt','Fzd3':'Wnt','Fzd8':'Wnt',
    'Lrp5':'Wnt','Lrp6':'Wnt','Ryk':'Wnt','Ror2':'Wnt',
    'Egfr':'EGF','Erbb2':'EGF',
    'Notch1':'Notch','Notch3':'Notch','Notch4':'Notch',
    'Tyro3':'TAM','Mertk':'TAM',
    'Nrp2':'Sema/VEGF','Plxna4':'Sema','Plxnb2':'Sema',
    'Epha3':'Ephrin','Epha7':'Ephrin','Ephb2':'Ephrin','Ephb3':'Ephrin','Ephb6':'Ephrin',
}
mouse_pairs['pathway'] = mouse_pairs['receptor'].map(pathway_map).fillna('Other')

print(f"Mouse mismatch pairs: {len(mouse_pairs)}")
print(f"\nTop 20:")
cols = ['receptor','ligand','r_score','r_logfc','l_score','l_logfc',
        'starvation','n_db','composite','pathway']
print(mouse_pairs[cols].head(20).to_string())

# ============================================================================
# 3. Cross-species conservation (using scores)
# ============================================================================
print("\n" + "=" * 80)
print("Step 3: Cross-species conservation")
print("=" * 80)

h_score = h_deg.set_index('gene')['scores'].to_dict()
h_padj = h_deg.set_index('gene')['pvals_adj'].to_dict()
h_logfc_map = h_deg.set_index('gene')['logfoldchanges'].to_dict()

results = []
# Deduplicate for cross-species iteration
lr_hc_unique = lr_hc[['receptor_human', 'ligand_human', 'receptor', 'ligand', 'n_db']].drop_duplicates(
    subset=['receptor_human', 'ligand_human'])

for _, lr_row in lr_hc_unique.iterrows():
    h_receptor = lr_row['receptor_human']
    h_ligand = lr_row['ligand_human']
    m_receptor = lr_row['receptor']
    m_ligand = lr_row['ligand']

    # Mouse
    m_r_s = m_score.get(m_receptor)
    m_r_p = m_padj.get(m_receptor)
    m_l_s = m_score.get(m_ligand)
    m_l_p = m_padj.get(m_ligand)
    m_r_lfc = m_logfc.get(m_receptor)
    m_l_lfc = m_logfc.get(m_ligand)

    # Human
    h_r_s = h_score.get(h_receptor)
    h_r_p = h_padj.get(h_receptor)
    h_l_s = h_score.get(h_ligand)
    h_l_p = h_padj.get(h_ligand)
    h_r_lfc = h_logfc_map.get(h_receptor)
    h_l_lfc = h_logfc_map.get(h_ligand)

    def check(r_s, r_p, l_s, l_p, strict=True):
        if r_s is None or l_s is None:
            return 'no_data'
        r_up = r_s > 0 and (r_p < 0.05 if strict else True)
        l_down = l_s < 0 and (l_p < 0.05 if strict else True)
        if r_up and l_down:
            return 'yes'
        elif r_s > 0 and l_s < 0:
            return 'trend'
        else:
            return 'no'

    m_strict = check(m_r_s, m_r_p, m_l_s, m_l_p, strict=True)
    m_relaxed = check(m_r_s, m_r_p, m_l_s, m_l_p, strict=False)
    h_strict = check(h_r_s, h_r_p, h_l_s, h_l_p, strict=True)
    h_relaxed = check(h_r_s, h_r_p, h_l_s, h_l_p, strict=False)

    # Conservation category
    if m_strict == 'yes' and h_strict == 'yes':
        conservation = 'CONSERVED_strict'
    elif m_relaxed == 'yes' and h_relaxed == 'yes':
        conservation = 'CONSERVED_relaxed'
    elif (m_strict == 'yes' or m_relaxed == 'yes') and h_relaxed in ('yes', 'trend'):
        conservation = 'CONSERVED_trend'
    elif m_strict == 'yes' or m_relaxed == 'yes':
        conservation = 'mouse_only'
    elif h_strict == 'yes' or h_relaxed == 'yes':
        conservation = 'human_only'
    else:
        conservation = 'neither'

    # Avg mismatch using scores (not logFC)
    m_mismatch = (min(m_r_s/200, 1) + min(abs(m_l_s)/200, 1)) / 2 if m_r_s and m_l_s and m_r_s > 0 and m_l_s < 0 else 0
    h_mismatch = (min(h_r_s/200, 1) + min(abs(h_l_s)/200, 1)) / 2 if h_r_s and h_l_s and h_r_s > 0 and h_l_s < 0 else 0

    results.append({
        'receptor': h_receptor, 'ligand': h_ligand, 'n_db': lr_row['n_db'],
        'm_r_score': round(m_r_s, 2) if m_r_s else None,
        'm_l_score': round(m_l_s, 2) if m_l_s else None,
        'm_r_lfc': round(m_r_lfc, 2) if m_r_lfc else None,
        'm_l_lfc': round(m_l_lfc, 2) if m_l_lfc else None,
        'h_r_score': round(h_r_s, 2) if h_r_s else None,
        'h_l_score': round(h_l_s, 2) if h_l_s else None,
        'h_r_lfc': round(h_r_lfc, 4) if h_r_lfc else None,
        'h_l_lfc': round(h_l_lfc, 4) if h_l_lfc else None,
        'm_strict': m_strict, 'h_strict': h_strict,
        'conservation': conservation,
        'm_mismatch_score': round(m_mismatch, 4),
        'h_mismatch_score': round(h_mismatch, 4),
        'avg_mismatch_score': round((m_mismatch + h_mismatch) / 2, 4),
    })

cross_df = pd.DataFrame(results)

print("Conservation summary:")
print(cross_df['conservation'].value_counts().to_string())

conserved = cross_df[cross_df['conservation'].str.startswith('CONSERVED')].copy()
conserved = conserved.sort_values('avg_mismatch_score', ascending=False).reset_index(drop=True)
print(f"\nConserved pairs: {len(conserved)}")

# ============================================================================
# 4. Final prioritization (using scores)
# ============================================================================
print("\n" + "=" * 80)
print("Step 4: Final prioritization")
print("=" * 80)

pairs = mouse_pairs.copy()
pairs['h_receptor'] = pairs['receptor'].str.upper()
pairs['h_ligand'] = pairs['ligand'].str.upper()

# Conservation
cons_map = {}
for _, row in cross_df.iterrows():
    cons_map[(row['receptor'], row['ligand'])] = row['conservation']

cons_scores = {
    'CONSERVED_strict': 1.0, 'CONSERVED_relaxed': 0.7,
    'CONSERVED_trend': 0.5, 'mouse_only': 0.2,
    'human_only': 0.1, 'neither': 0.0, 'not_checked': 0.1,
}
pairs['conservation'] = pairs.apply(
    lambda r: cons_map.get((r['h_receptor'], r['h_ligand']), 'not_checked'), axis=1)
pairs['conservation_score'] = pairs['conservation'].map(cons_scores)

# Druggability (corrected: UniProt + DGIdb)
pairs['lig_deliver'] = pairs['h_ligand'].map(
    lambda g: deliver_data.loc[g, 'deliver_score'] if g in deliver_data.index else 0.3)
pairs['lig_dgidb'] = pairs['h_ligand'].map(
    lambda g: gene_scores.loc[g, 'druggability_auto'] if g in gene_scores.index else 0.1)
pairs['rec_dgidb'] = pairs['h_receptor'].map(
    lambda g: gene_scores.loc[g, 'druggability_auto'] if g in gene_scores.index else 0.1)
pairs['druggability'] = 0.3 * pairs['rec_dgidb'] + 0.7 * (0.6 * pairs['lig_deliver'] + 0.4 * pairs['lig_dgidb'])

# Literature
pairs['lit_receptor'] = pairs['h_receptor'].map(
    lambda g: gene_scores.loc[g, 'literature_auto'] if g in gene_scores.index else 0.0)
pairs['lit_ligand'] = pairs['h_ligand'].map(
    lambda g: gene_scores.loc[g, 'literature_auto'] if g in gene_scores.index else 0.0)
pairs['literature'] = (pairs['lit_receptor'] + pairs['lit_ligand']) / 2

# Mismatch norm
max_m = pairs['composite'].max()
pairs['mismatch_norm'] = pairs['composite'] / max_m if max_m > 0 else 0

# Final score
W_M, W_C, W_D, W_L = 0.30, 0.30, 0.20, 0.20
pairs['priority_score'] = (
    W_M * pairs['mismatch_norm'] +
    W_C * pairs['conservation_score'] +
    W_D * pairs['druggability'] +
    W_L * pairs['literature']
)

pairs = pairs.sort_values('priority_score', ascending=False).reset_index(drop=True)
pairs.index += 1
pairs.index.name = 'rank'

# ============================================================================
# 5. Display results
# ============================================================================
print("\nFINAL RANKING (Wilcoxon scores)")

print(f"\n{'Rank':<5} {'Pair':<22} {'TOTAL':>6} = "
      f"{'Mismatch':>8} + {'Conserv':>7} + {'Drug':>6} + {'Lit':>6}  "
      f"{'R_score':>8} {'L_score':>8} {'Pathway'}")
print("-" * 110)
for i in range(min(20, len(pairs))):
    r = pairs.iloc[i]
    pair = f"{r['receptor']}/{r['ligand']}"
    print(f" {i+1:<4} {pair:<22} {r['priority_score']:.3f} = "
          f"{W_M*r['mismatch_norm']:>8.3f} + {W_C*r['conservation_score']:>7.3f} + "
          f"{W_D*r['druggability']:>6.3f} + {W_L*r['literature']:>6.3f}  "
          f"{r['r_score']:>8.1f} {r['l_score']:>8.1f} {r.get('pathway','')}")

# Key pairs comparison
print("\n" + "=" * 80)
print("Key pairs: logFC-based vs score-based ranking")
print("=" * 80)

old = pd.read_csv(MISMATCH_DIR / "therapeutic_targets_corrected.csv", index_col=0)

key = [('Fgfr2','Fgf10'),('Fgfr2','Fgf7'),('Acvr1','Bmp6'),('Bmpr2','Bmp6'),
       ('Tyro3','Gas6'),('Insr','Nampt'),('Unc5b','Ntn1'),('Epha7','Efna2'),
       ('Bmpr2','Bmp2'),('Notch1','Psen1')]

print(f"\n{'Pair':<22} {'logFC rank':>12} {'Score rank':>12} {'Change':>8}")
print("-" * 58)
for rec, lig in key:
    o = old[(old['receptor']==rec) & (old['ligand']==lig)]
    n = pairs[(pairs['receptor']==rec) & (pairs['ligand']==lig)]
    o_r = o.index[0] if len(o) > 0 else '-'
    n_r = n.index[0] if len(n) > 0 else '-'
    if isinstance(o_r, int) and isinstance(n_r, int):
        change = o_r - n_r
        arrow = f"{'↑' if change > 0 else '↓'}{abs(change)}" if change != 0 else '='
    else:
        arrow = ''
    print(f"  {rec}/{lig:<14} {str(o_r):>12} {str(n_r):>12} {arrow:>8}")

# ============================================================================
# 6. Save
# ============================================================================
print("\n" + "=" * 80)
print("Saving")
print("=" * 80)

mouse_pairs.to_csv(OUTPUT_DIR / "mouse_lr_mismatch_scores.csv")
cross_df.to_csv(OUTPUT_DIR / "cross_species_lr_mismatch_scores.csv", index=False)
conserved.to_csv(OUTPUT_DIR / "cross_species_conserved_scores.csv", index=False)
pairs.to_csv(OUTPUT_DIR / "therapeutic_targets_scores.csv")

print(f"Saved: mouse_lr_mismatch_scores.csv ({len(mouse_pairs)} pairs)")
print(f"Saved: cross_species_lr_mismatch_scores.csv ({len(cross_df)} pairs)")
print(f"Saved: cross_species_conserved_scores.csv ({len(conserved)} pairs)")
print(f"Saved: therapeutic_targets_scores.csv ({len(pairs)} pairs)")
