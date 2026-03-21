#!/usr/bin/env python3
"""
Phase 6: Therapeutic Target Prioritization

Scoring (Geneformer skipped, weights redistributed):
  1. Mismatch Score (composite from refined analysis)  — 30%
  2. Cross-species Conservation                         — 30%
  3. Druggability (recombinant ligand available)        — 20%
  4. Literature Support (cardiac/epicardial relevance)  — 20%

Input: mouse refined mismatch + cross-species conservation data
Output: Final ranked therapeutic target list
"""

import pandas as pd
import numpy as np
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent.parent.parent
MISMATCH_DIR = PROJECT_DIR / "results" / "mismatch"
OUTPUT_DIR = MISMATCH_DIR

# ============================================================================
# 1. Load data
# ============================================================================
print("=" * 80)
print("Step 1: Loading mismatch and conservation data")
print("=" * 80)

mouse_refined = pd.read_csv(MISMATCH_DIR / "mouse_lr_mismatch_refined.csv")
cross_species = pd.read_csv(MISMATCH_DIR / "cross_species_conserved.csv")
cross_all = pd.read_csv(MISMATCH_DIR / "cross_species_lr_mismatch.csv")

print(f"Mouse refined pairs: {len(mouse_refined)}")
print(f"Cross-species conserved: {len(cross_species)}")
print(f"Cross-species all: {len(cross_all)}")

# ============================================================================
# 2. Build unified pair table from mouse refined results
# ============================================================================
print("\n" + "=" * 80)
print("Step 2: Building unified pair table")
print("=" * 80)

# Start with mouse refined pairs (127 pairs, already high-quality)
pairs = mouse_refined[['receptor', 'ligand', 'r_logfc', 'l_logfc', 'l_padj',
                        'n_db', 'starvation', 'composite', 'pathway']].copy()
pairs = pairs.rename(columns={'composite': 'mismatch_composite'})

# Convert mouse gene names to human for matching cross-species data
pairs['h_receptor'] = pairs['receptor'].str.upper()
pairs['h_ligand'] = pairs['ligand'].str.upper()

# ============================================================================
# 3. Add conservation scores
# ============================================================================
print("\n" + "=" * 80)
print("Step 3: Adding conservation scores")
print("=" * 80)

# Build conservation lookup from cross-species data
conservation_map = {}
for _, row in cross_all.iterrows():
    key = (row['receptor'], row['ligand'])
    conservation_map[key] = row['conservation']

def get_conservation(row):
    key = (row['h_receptor'], row['h_ligand'])
    return conservation_map.get(key, 'not_checked')

pairs['conservation'] = pairs.apply(get_conservation, axis=1)

# Conservation score
conservation_scores = {
    'CONSERVED_strict': 1.0,
    'CONSERVED_relaxed': 0.7,
    'CONSERVED_trend': 0.5,
    'mouse_only': 0.2,
    'human_only': 0.1,
    'neither': 0.0,
    'not_checked': 0.1,  # not in HC database
}
pairs['conservation_score'] = pairs['conservation'].map(conservation_scores)

print(pairs['conservation'].value_counts().to_string())

# ============================================================================
# 4. Add druggability scores
# ============================================================================
print("\n" + "=" * 80)
print("Step 4: Adding druggability scores")
print("=" * 80)

# Druggability: is the ligand available as a recombinant protein or approved drug?
# Based on known recombinant proteins and drug availability
druggability = {
    # FGF family — recombinant FGFs widely available
    'Fgf1': 1.0, 'Fgf2': 1.0, 'Fgf7': 1.0, 'Fgf10': 1.0, 'Fgf9': 0.8,
    # BMP family — recombinant BMPs available (BMP2/7 FDA-approved for bone)
    'Bmp2': 1.0, 'Bmp4': 0.9, 'Bmp6': 0.8, 'Bmp7': 1.0,
    'Inhbb': 0.5,  # Activin/inhibin — research grade
    # Wnt — difficult to produce, some small molecule agonists
    'Wnt3': 0.3, 'Wnt5a': 0.4, 'Wnt9a': 0.3, 'Wnt5b': 0.3,
    'Wnt4': 0.3, 'Rspo1': 0.6, 'Rspo2': 0.5, 'Dkk1': 0.7,
    # TGF-beta — recombinant available
    'Tgfb1': 0.9,
    # Growth factors — recombinant available
    'Gas6': 0.7, 'Ntn1': 0.6, 'Pgf': 0.7,
    # Semaphorins — research grade
    'Sema4c': 0.4, 'Sema6a': 0.4, 'Sema3c': 0.4, 'Sema3f': 0.4,
    # Ephrins — Fc-fused proteins available
    'Efna1': 0.6, 'Efna2': 0.6, 'Efnb2': 0.6,
    'Afdn': 0.1,  # intracellular, not druggable as ligand
    # Interleukins — recombinant widely available
    'Il2': 1.0, 'Il4': 0.9, 'Il18': 0.8, 'Il1b': 0.9, 'Il1rn': 1.0,
    'Il23a': 0.7,
    # Neuropeptides
    'Penk': 0.5, 'Tslp': 0.7, 'Edn2': 0.6, 'Edn3': 0.6,
    # ECM proteins
    'Vtn': 0.5, 'Ccn1': 0.4, 'Ltbp1': 0.3,
    # Small molecules can target receptors directly
    'Nampt': 0.6,  # NAD biosynthesis enzyme, inhibitors exist
}

pairs['druggability'] = pairs['ligand'].map(druggability).fillna(0.2)
print(f"Ligands with known druggability: {pairs['druggability'].gt(0.2).sum()}")

# ============================================================================
# 5. Add literature support scores
# ============================================================================
print("\n" + "=" * 80)
print("Step 5: Adding literature support scores")
print("=" * 80)

# Literature support for receptor-ligand axis in cardiac regeneration/epicardial biology
# Scored 0-1 based on existing evidence
literature = {
    # FGF — extensive evidence in epicardial activation and cardiac repair
    ('Fgfr2', 'Fgf10'): 1.0,  # Cheng lab validated, multiple papers
    ('Fgfr2', 'Fgf7'): 0.7,   # KGF, known FGFR2b ligand, some cardiac evidence
    ('Fgfr2', 'Fgf1'): 0.6,   # aFGF, general cardiac evidence
    ('Fgfr2', 'Fgf2'): 0.6,   # bFGF, extensive cardiac literature
    ('Fgfrl1', 'Fgf1'): 0.3,  # decoy receptor, less studied
    ('Fgfrl1', 'Fgf2'): 0.3,
    # BMP — substantial evidence in cardiac development and repair
    ('Bmpr2', 'Bmp4'): 0.7,   # cardiac development
    ('Bmpr2', 'Bmp6'): 0.5,
    ('Bmpr2', 'Bmp2'): 0.8,   # BMP2 well-studied in heart
    ('Bmpr1a', 'Bmp4'): 0.7,
    ('Bmpr1a', 'Bmp6'): 0.5,
    ('Bmpr1a', 'Bmp2'): 0.7,
    ('Acvr1', 'Bmp6'): 0.5,
    ('Acvr1', 'Inhbb'): 0.4,
    ('Acvr1', 'Bmp2'): 0.6,
    # Wnt — important in cardiac development, some epicardial evidence
    ('Lrp6', 'Rspo1'): 0.5,
    ('Lrp6', 'Dkk1'): 0.5,
    ('Fzd1', 'Wnt3'): 0.3,
    ('Fzd8', 'Wnt3'): 0.3,
    # Notch — cardiac development and regeneration
    ('Notch1', 'Mfng'): 0.4,
    ('Notch3', 'Psen1'): 0.4,
    # Netrin/Guidance — some cardiac innervation evidence
    ('Unc5b', 'Ntn1'): 0.4,
    ('Neo1', 'Ntn1'): 0.3,
    ('Neo1', 'Rgma'): 0.2,
    ('Neo1', 'Rgmb'): 0.2,
    # Ephrin — EMT and cell migration, some cardiac evidence
    ('Epha7', 'Efna1'): 0.3,
    ('Epha7', 'Efna2'): 0.3,
    ('Epha3', 'Efna2'): 0.3,
    ('Ephb6', 'Afdn'): 0.2,
    # TAM receptors — efferocytosis post-MI
    ('Tyro3', 'Gas6'): 0.5,
    ('Mertk', 'Gas6'): 0.6,  # MERTK well-studied in cardiac macrophages
    # TGF-beta — extensive cardiac fibrosis literature
    ('Tgfbr1', 'Tgfb1'): 0.8,
    # EGF
    ('Egfr', 'Efemp1'): 0.2,
    # Integrin
    ('Itgb8', 'Vtn'): 0.2,
    ('Itgb8', 'Tgfb1'): 0.4,
    # VEGF/NRP
    ('Nrp2', 'Pgf'): 0.4,
    ('Nrp2', 'Sema3c'): 0.3,
    # Insulin receptor
    ('Insr', 'Nampt'): 0.3,
    # Endothelin
    ('Ednra', 'Edn3'): 0.4,
    # Adenosine
    ('Adora2b', 'Ntn1'): 0.2,
}

pairs['literature'] = pairs.apply(
    lambda r: literature.get((r['receptor'], r['ligand']), 0.1), axis=1)
print(f"Pairs with literature annotation: {pairs['literature'].gt(0.1).sum()}")

# ============================================================================
# 6. Calculate final priority score
# ============================================================================
print("\n" + "=" * 80)
print("Step 6: Calculating final priority score")
print("=" * 80)

# Normalize mismatch composite to [0, 1]
max_mismatch = pairs['mismatch_composite'].max()
pairs['mismatch_norm'] = pairs['mismatch_composite'] / max_mismatch if max_mismatch > 0 else 0

# Final score (weights redistributed without Geneformer)
W_MISMATCH = 0.30
W_CONSERVATION = 0.30
W_DRUGGABILITY = 0.20
W_LITERATURE = 0.20

pairs['priority_score'] = (
    W_MISMATCH * pairs['mismatch_norm'] +
    W_CONSERVATION * pairs['conservation_score'] +
    W_DRUGGABILITY * pairs['druggability'] +
    W_LITERATURE * pairs['literature']
)

pairs = pairs.sort_values('priority_score', ascending=False).reset_index(drop=True)
pairs.index += 1
pairs.index.name = 'rank'

# ============================================================================
# 7. Display results
# ============================================================================
print("\n" + "=" * 80)
print("FINAL THERAPEUTIC TARGET RANKING")
print(f"Scoring: Mismatch {W_MISMATCH:.0%} + Conservation {W_CONSERVATION:.0%} "
      f"+ Druggability {W_DRUGGABILITY:.0%} + Literature {W_LITERATURE:.0%}")
print("=" * 80)

display_cols = ['receptor', 'ligand', 'priority_score',
                'mismatch_norm', 'conservation', 'conservation_score',
                'druggability', 'literature', 'n_db', 'pathway']

print("\nTop 30:")
print(pairs[display_cols].head(30).to_string())

# Score breakdown for top 15
print("\n" + "=" * 80)
print("Score breakdown (top 15)")
print("=" * 80)
print(f"\n{'Rank':<5} {'Pair':<22} {'TOTAL':>6} = "
      f"{'Mismatch':>8}({W_MISMATCH:.0%}) + {'Conserv':>7}({W_CONSERVATION:.0%}) + "
      f"{'Drug':>5}({W_DRUGGABILITY:.0%}) + {'Lit':>5}({W_LITERATURE:.0%})  {'Pathway':<15}")
print("-" * 105)
for i in range(min(15, len(pairs))):
    row = pairs.iloc[i]
    rank = i + 1
    pair_name = f"{row['receptor']}/{row['ligand']}"
    print(f" {rank:<4} {pair_name:<22} {row['priority_score']:>5.3f} = "
          f"{W_MISMATCH*row['mismatch_norm']:>8.3f}      + {W_CONSERVATION*row['conservation_score']:>7.3f}       + "
          f"{W_DRUGGABILITY*row['druggability']:>5.3f}     + {W_LITERATURE*row['literature']:>5.3f}      {row.get('pathway',''):>15}")

# ============================================================================
# 8. Validation: Where is FGF10/FGFR2?
# ============================================================================
print("\n" + "=" * 80)
print("Validation: FGF10/FGFR2 and instruction 3 expected targets")
print("=" * 80)

expected = [
    ('Fgfr2', 'Fgf10', 'FGF10 delivery'),
    ('Bmpr2', 'Bmp2', 'BMP2 delivery (FDA-approved)'),
    ('Bmpr2', 'Bmp4', 'BMP4 delivery'),
    ('Acvr1', 'Bmp6', 'BMP6 delivery'),
    ('Fgfr2', 'Fgf7', 'FGF7/KGF delivery'),
    ('Fgfr2', 'Fgf1', 'FGF1 delivery'),
    ('Tgfbr1', 'Tgfb1', 'TGF-β1 modulation'),
]

for receptor, ligand, action in expected:
    m = pairs[(pairs['receptor']==receptor) & (pairs['ligand']==ligand)]
    if len(m) > 0:
        rank = m.index[0]
        row = m.iloc[0]
        print(f"  {receptor}/{ligand:<10} Rank {rank:>3}/{len(pairs)}  "
              f"score={row['priority_score']:.3f}  conserv={row['conservation']:<20}  → {action}")
    else:
        print(f"  {receptor}/{ligand:<10} not in refined set")

# ============================================================================
# 9. Suggested therapeutic candidates (actionable)
# ============================================================================
print("\n" + "=" * 80)
print("ACTIONABLE THERAPEUTIC CANDIDATES")
print("(Priority score > 0.3, druggability > 0.5)")
print("=" * 80)

actionable = pairs[(pairs['priority_score'] > 0.3) & (pairs['druggability'] >= 0.5)]
if len(actionable) > 0:
    cols = ['receptor', 'ligand', 'priority_score', 'conservation',
            'druggability', 'literature', 'pathway']
    print(actionable[cols].to_string())
else:
    print("No pairs meet both criteria. Showing top 10 by priority:")
    print(pairs[['receptor', 'ligand', 'priority_score', 'conservation',
                 'druggability', 'literature', 'pathway']].head(10).to_string())

# ============================================================================
# 10. Save
# ============================================================================
print("\n" + "=" * 80)
print("Saving")
print("=" * 80)

pairs.to_csv(OUTPUT_DIR / "therapeutic_targets_prioritized.csv")
print(f"Saved: therapeutic_targets_prioritized.csv ({len(pairs)} pairs)")

if len(actionable) > 0:
    actionable.to_csv(OUTPUT_DIR / "therapeutic_targets_actionable.csv")
    print(f"Saved: therapeutic_targets_actionable.csv ({len(actionable)} pairs)")
