#!/usr/bin/env python3
"""
Fix druggability bias by adding UniProt subcellular location.

DGIdb captures traditional drug-target interactions (small molecules, antibodies)
but misses recombinant protein delivery — our primary therapeutic strategy.

Fix: Query UniProt for subcellular location of each ligand.
  - Secreted → high deliverability as recombinant protein
  - Cell membrane → medium (Fc-fusion possible)
  - Intracellular → low

New druggability = weighted combination of:
  - DGIdb score (traditional drugs, targets receptor)
  - UniProt deliverability (recombinant protein, targets ligand)
"""

import pandas as pd
import numpy as np
import json
import urllib.request
import urllib.parse
import time
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent.parent.parent
MISMATCH_DIR = PROJECT_DIR / "results" / "mismatch"

# ============================================================================
# 1. Load existing automated scores
# ============================================================================
print("=" * 80)
print("Step 1: Loading existing scores")
print("=" * 80)

gene_scores = pd.read_csv(MISMATCH_DIR / "gene_scores_automated.csv", index_col=0)
print(f"Genes with existing scores: {len(gene_scores)}")

# ============================================================================
# 2. Query UniProt for subcellular location
# ============================================================================
print("\n" + "=" * 80)
print("Step 2: Querying UniProt for subcellular location")
print("=" * 80)

def query_uniprot_location(gene):
    """Query UniProt for subcellular location of a human gene."""
    params = urllib.parse.urlencode({
        'query': f'(gene_exact:{gene}) AND (organism_id:9606) AND (reviewed:true)',
        'fields': 'accession,gene_names,cc_subcellular_location',
        'format': 'json',
        'size': '1'
    })
    url = f'https://rest.uniprot.org/uniprotkb/search?{params}'

    req = urllib.request.Request(url, headers={
        'User-Agent': 'EpicardialAnalysis/1.0',
        'Accept': 'application/json'
    })
    resp = urllib.request.urlopen(req, timeout=10)
    data = json.loads(resp.read())

    if not data['results']:
        return []

    locations = []
    for comment in data['results'][0].get('comments', []):
        if comment.get('commentType') == 'SUBCELLULAR LOCATION':
            for sub in comment.get('subcellularLocations', []):
                loc = sub.get('location', {}).get('value', '')
                if loc:
                    locations.append(loc)
    return locations


all_genes = gene_scores.index.tolist()
uniprot_results = {}

for i, gene in enumerate(all_genes):
    try:
        locs = query_uniprot_location(gene)
        uniprot_results[gene] = locs
    except Exception as e:
        uniprot_results[gene] = []

    if (i + 1) % 3 == 0:
        time.sleep(0.35)
    if (i + 1) % 30 == 0:
        print(f"  Queried {i+1}/{len(all_genes)} genes...")

print(f"\nUniProt results: {len(uniprot_results)} genes queried")

# ============================================================================
# 3. Classify protein deliverability
# ============================================================================
print("\n" + "=" * 80)
print("Step 3: Classifying protein deliverability")
print("=" * 80)

def classify_deliverability(locations):
    """
    Score deliverability based on subcellular location.
    Secreted proteins are ideal for recombinant delivery.
    """
    if not locations:
        return 'unknown', 0.3

    loc_str = ' '.join(locations).lower()

    if 'secreted' in loc_str:
        return 'secreted', 1.0
    elif 'cell membrane' in loc_str or 'membrane' in loc_str:
        # Membrane proteins: extracellular domain can be delivered as Fc-fusion
        # or the protein can be a membrane-bound ligand (e.g., ephrins)
        if any(x in loc_str for x in ['single-pass type i', 'gpi-anchor', 'lipid-anchor']):
            return 'membrane_ectodomain', 0.6
        else:
            return 'membrane', 0.5
    elif 'cytoplasm' in loc_str or 'nucleus' in loc_str:
        return 'intracellular', 0.1
    elif 'extracellular' in loc_str:
        return 'extracellular', 0.8
    else:
        return 'other', 0.3

deliver_data = {}
for gene, locs in uniprot_results.items():
    category, score = classify_deliverability(locs)
    deliver_data[gene] = {
        'uniprot_locations': ' | '.join(locs) if locs else 'N/A',
        'deliver_category': category,
        'deliver_score': score,
    }

deliver_df = pd.DataFrame(deliver_data).T
deliver_df.index.name = 'gene'

# Summary
print("Deliverability category distribution:")
print(deliver_df['deliver_category'].value_counts().to_string())
print()

# Show key genes
for gene in ['FGF10', 'FGF7', 'BMP6', 'BMP4', 'GAS6', 'EFNB2',
             'HRAS', 'PSEN1', 'NTN1', 'EDN3', 'NAMPT', 'IL2']:
    if gene in deliver_data:
        d = deliver_data[gene]
        print(f"  {gene:<10} {d['deliver_category']:<22} score={d['deliver_score']:.1f}  "
              f"locs={d['uniprot_locations'][:60]}")

# ============================================================================
# 4. Compute corrected druggability
# ============================================================================
print("\n" + "=" * 80)
print("Step 4: Computing corrected druggability scores")
print("=" * 80)

# Merge into gene_scores
gene_scores['deliver_category'] = gene_scores.index.map(
    lambda g: deliver_data.get(g, {}).get('deliver_category', 'unknown'))
gene_scores['deliver_score'] = gene_scores.index.map(
    lambda g: deliver_data.get(g, {}).get('deliver_score', 0.3))
gene_scores['uniprot_locations'] = gene_scores.index.map(
    lambda g: deliver_data.get(g, {}).get('uniprot_locations', 'N/A'))

# Corrected druggability:
# For the "deliver ligand" strategy:
#   - Ligand deliverability (UniProt) matters most (weight 0.6)
#   - Traditional drug availability (DGIdb) as bonus (weight 0.4)
# This is applied at the pair level, not gene level

# ============================================================================
# 5. Re-run prioritization
# ============================================================================
print("\n" + "=" * 80)
print("Step 5: Re-running prioritization with corrected druggability")
print("=" * 80)

mouse_refined = pd.read_csv(MISMATCH_DIR / "mouse_lr_mismatch_refined.csv")
cross_all = pd.read_csv(MISMATCH_DIR / "cross_species_lr_mismatch.csv")

pairs = mouse_refined[['receptor', 'ligand', 'r_logfc', 'l_logfc',
                        'n_db', 'starvation', 'composite', 'pathway']].copy()
pairs = pairs.rename(columns={'composite': 'mismatch_composite'})
pairs['h_receptor'] = pairs['receptor'].str.upper()
pairs['h_ligand'] = pairs['ligand'].str.upper()

# Conservation
conservation_map = {}
for _, row in cross_all.iterrows():
    conservation_map[(row['receptor'], row['ligand'])] = row['conservation']

conservation_scores = {
    'CONSERVED_strict': 1.0, 'CONSERVED_relaxed': 0.7,
    'CONSERVED_trend': 0.5, 'mouse_only': 0.2,
    'human_only': 0.1, 'neither': 0.0, 'not_checked': 0.1,
}
pairs['conservation'] = pairs.apply(
    lambda r: conservation_map.get((r['h_receptor'], r['h_ligand']), 'not_checked'), axis=1)
pairs['conservation_score'] = pairs['conservation'].map(conservation_scores)

# Corrected druggability:
# Ligand: 60% deliverability (UniProt) + 40% DGIdb
# Then weight: 30% receptor DGIdb + 70% ligand combined
pairs['lig_deliver'] = pairs['h_ligand'].map(
    lambda g: deliver_data.get(g, {}).get('deliver_score', 0.3))
pairs['lig_dgidb'] = pairs['h_ligand'].map(
    lambda g: gene_scores.loc[g, 'druggability_auto'] if g in gene_scores.index else 0.1)
pairs['rec_dgidb'] = pairs['h_receptor'].map(
    lambda g: gene_scores.loc[g, 'druggability_auto'] if g in gene_scores.index else 0.1)

pairs['lig_drug_combined'] = 0.6 * pairs['lig_deliver'] + 0.4 * pairs['lig_dgidb']
pairs['druggability_corrected'] = 0.3 * pairs['rec_dgidb'] + 0.7 * pairs['lig_drug_combined']

# Literature (same as before)
pairs['lit_receptor'] = pairs['h_receptor'].map(
    lambda g: gene_scores.loc[g, 'literature_auto'] if g in gene_scores.index else 0.0)
pairs['lit_ligand'] = pairs['h_ligand'].map(
    lambda g: gene_scores.loc[g, 'literature_auto'] if g in gene_scores.index else 0.0)
pairs['literature_auto'] = (pairs['lit_receptor'] + pairs['lit_ligand']) / 2

# Mismatch norm
max_m = pairs['mismatch_composite'].max()
pairs['mismatch_norm'] = pairs['mismatch_composite'] / max_m if max_m > 0 else 0

# Final score
W_M, W_C, W_D, W_L = 0.30, 0.30, 0.20, 0.20

pairs['priority_score'] = (
    W_M * pairs['mismatch_norm'] +
    W_C * pairs['conservation_score'] +
    W_D * pairs['druggability_corrected'] +
    W_L * pairs['literature_auto']
)

pairs = pairs.sort_values('priority_score', ascending=False).reset_index(drop=True)
pairs.index += 1
pairs.index.name = 'rank'

# ============================================================================
# 6. Display
# ============================================================================
cols = ['receptor', 'ligand', 'priority_score', 'mismatch_norm',
        'conservation_score', 'druggability_corrected', 'lig_deliver',
        'literature_auto', 'n_db', 'pathway']

print("\nTop 30:")
print(pairs[cols].head(30).to_string())

# Score breakdown
print(f"\n{'Rank':<5} {'Pair':<22} {'TOTAL':>6} = "
      f"{'Mismatch':>8} + {'Conserv':>7} + {'Drug_fix':>8} + {'Lit':>6}  "
      f"{'Deliver':<10} {'Pathway'}")
print("-" * 105)
for i in range(min(20, len(pairs))):
    r = pairs.iloc[i]
    pair = f"{r['receptor']}/{r['ligand']}"
    cat = deliver_data.get(r['h_ligand'], {}).get('deliver_category', '?')
    print(f" {i+1:<4} {pair:<22} {r['priority_score']:.3f} = "
          f"{W_M*r['mismatch_norm']:>8.3f} + {W_C*r['conservation_score']:>7.3f} + "
          f"{W_D*r['druggability_corrected']:>8.3f} + {W_L*r['literature_auto']:>6.3f}  "
          f"{cat:<10} {r.get('pathway','')}")

# Comparison: manual vs DGIdb-only vs corrected
print("\n" + "=" * 80)
print("Comparison across scoring methods")
print("=" * 80)

# Load previous results for comparison
manual = pd.read_csv(MISMATCH_DIR / "therapeutic_targets_prioritized.csv", index_col=0)
auto = pd.read_csv(MISMATCH_DIR / "therapeutic_targets_automated.csv", index_col=0)

key_pairs = [('Fgfr2','Fgf10'),('Fgfr2','Fgf7'),('Acvr1','Bmp6'),
             ('Bmpr2','Bmp6'),('Tyro3','Gas6'),('Insr','Nampt'),
             ('Unc5b','Ntn1'),('Bmpr2','Bmp2'),('Epha7','Efna2')]

print(f"\n{'Pair':<22} {'Manual':>8} {'DGIdb-only':>12} {'Corrected':>12}")
print("-" * 58)
for rec, lig in key_pairs:
    m_rank = manual[(manual['receptor']==rec) & (manual['ligand']==lig)]
    a_rank = auto[(auto['receptor']==rec) & (auto['ligand']==lig)]
    c_rank = pairs[(pairs['receptor']==rec) & (pairs['ligand']==lig)]

    m_r = m_rank.index[0] if len(m_rank) > 0 else '-'
    a_r = a_rank.index[0] if len(a_rank) > 0 else '-'
    c_r = c_rank.index[0] if len(c_rank) > 0 else '-'
    print(f"  {rec}/{lig:<14} {str(m_r):>8} {str(a_r):>12} {str(c_r):>12}")

# ============================================================================
# 7. Save
# ============================================================================
print("\n" + "=" * 80)
print("Saving")
print("=" * 80)

pairs.to_csv(MISMATCH_DIR / "therapeutic_targets_corrected.csv")
gene_scores.to_csv(MISMATCH_DIR / "gene_scores_corrected.csv")
deliver_df.to_csv(MISMATCH_DIR / "uniprot_deliverability.csv")

print(f"Saved: therapeutic_targets_corrected.csv ({len(pairs)} pairs)")
print(f"Saved: gene_scores_corrected.csv ({len(gene_scores)} genes)")
print(f"Saved: uniprot_deliverability.csv ({len(deliver_df)} genes)")
