#!/usr/bin/env python3
"""
Phase 1: Kuppe 2022 Data Analysis
Project: AI-Driven Epicardial Target Discovery in MI

Analyze myocardial infarction single-cell data to identify epicardial cells
"""

import anndata as ad
import numpy as np
import pandas as pd
import os
from pathlib import Path

# Paths
PROJECT_DIR = Path(__file__).parent.parent.parent
DATA_FILE = PROJECT_DIR / "data/raw/kuppe/54d24dbe-d39a-4844-bb21-07b5f4e173ad.h5ad"
OUTPUT_DIR = PROJECT_DIR / "data/processed"
RESULTS_DIR = PROJECT_DIR / "results/tables"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Epicardial markers (TCF21 included for reference but noted as low expression)
EPICARDIAL_MARKERS = {
    'WT1': 'ENSG00000184937',
    'TBX18': 'ENSG00000112837',
    'ALDH1A2': 'ENSG00000128918',
    'UPK3B': 'ENSG00000243566',
    'TCF21': 'ENSG00000118526'  # Often low in adult epicardium
}

def analyze_data():
    print("="*60)
    print("Kuppe 2022 MI Heart Data Analysis")
    print("="*60)

    # Load with backed mode for memory efficiency
    print("\n1. Loading data (backed mode)...")
    adata = ad.read_h5ad(DATA_FILE, backed='r')

    print(f"   Cells: {adata.n_obs:,}")
    print(f"   Genes: {adata.n_vars:,}")

    # Basic statistics
    print("\n2. Cell type distribution:")
    cell_counts = adata.obs['cell_type_original'].value_counts()
    for ct, count in cell_counts.items():
        pct = 100 * count / adata.n_obs
        print(f"   {ct}: {count:,} ({pct:.1f}%)")

    print("\n3. Disease state distribution:")
    disease_counts = adata.obs['disease'].value_counts()
    for d, count in disease_counts.items():
        pct = 100 * count / adata.n_obs
        print(f"   {d}: {count:,} ({pct:.1f}%)")

    print("\n4. Patient groups:")
    group_counts = adata.obs['patient_group'].value_counts()
    for g, count in group_counts.items():
        pct = 100 * count / adata.n_obs
        print(f"   {g}: {count:,} ({pct:.1f}%)")

    # Check epicardial markers
    print("\n5. Epicardial marker genes:")
    markers_found = {}
    for name, ensembl in EPICARDIAL_MARKERS.items():
        if ensembl in adata.var_names:
            print(f"   + {name}: {ensembl}")
            markers_found[name] = ensembl
        else:
            print(f"   - {name}: NOT FOUND")

    # Save metadata summary
    print("\n6. Saving metadata summary...")

    # Cell type summary
    summary_df = pd.DataFrame({
        'cell_type': cell_counts.index,
        'count': cell_counts.values,
        'percentage': (100 * cell_counts.values / adata.n_obs).round(2)
    })
    summary_df.to_csv(RESULTS_DIR / "kuppe_celltype_summary.csv", index=False)

    # Sample metadata
    sample_meta = adata.obs[['sample', 'donor_id', 'patient_group', 'disease',
                             'cell_type_original']].drop_duplicates(subset=['sample'])
    sample_meta.to_csv(RESULTS_DIR / "kuppe_sample_metadata.csv", index=False)

    # Full obs metadata (useful for downstream)
    print("   Saving cell metadata...")
    cell_meta = adata.obs.copy()
    cell_meta.to_csv(OUTPUT_DIR / "kuppe_cell_metadata.csv")

    # Gene metadata
    gene_meta = adata.var.copy()
    gene_meta['ensembl_id'] = adata.var_names
    gene_meta.to_csv(OUTPUT_DIR / "kuppe_gene_metadata.csv")

    print(f"\n   Saved to: {RESULTS_DIR}")
    print(f"   Saved to: {OUTPUT_DIR}")

    return adata, markers_found

def identify_epicardial_by_annotation(adata):
    """Identify epicardial cells by existing annotations"""
    print("\n" + "="*60)
    print("Identifying Epicardial Cells")
    print("="*60)

    # Look for epicardial-related annotations
    epi_keywords = ['epicard', 'meso', 'adipocyte']

    cell_types = adata.obs['cell_type'].unique()
    orig_types = adata.obs['cell_type_original'].unique()

    print("\n1. Searching for epicardial-related cell types...")

    epi_types = []
    for ct in set(cell_types) | set(orig_types):
        ct_lower = str(ct).lower()
        if any(kw in ct_lower for kw in epi_keywords):
            epi_types.append(ct)
            print(f"   Found: {ct}")

    # Count cells
    epi_mask = adata.obs['cell_type_original'].isin(epi_types)
    n_epi = epi_mask.sum()
    print(f"\n2. Cells with epicardial-related annotations: {n_epi:,} ({100*n_epi/adata.n_obs:.2f}%)")

    if n_epi > 0:
        print("\n   Breakdown:")
        for ct in epi_types:
            count = (adata.obs['cell_type_original'] == ct).sum()
            if count > 0:
                print(f"     {ct}: {count:,}")

    # Note: In this dataset, "Adipocyte" cells are actually epicardial fat adipocytes
    # True epicardial cells (mesothelial) might be mixed with fibroblasts

    print("\n3. Note: This dataset annotates 'Adipocyte' as epicardial fat cells.")
    print("   True epicardial mesothelial cells may be within 'Fibroblast' population.")
    print("   Will need marker-based identification for refined analysis.")

    return epi_mask

def main():
    print("\n" + "#"*60)
    print("# Phase 1: Kuppe Data Exploration")
    print("#"*60 + "\n")

    # Analyze data
    adata, markers = analyze_data()

    # Identify epicardial cells
    epi_mask = identify_epicardial_by_annotation(adata)

    # Summary
    print("\n" + "="*60)
    print("Summary")
    print("="*60)
    print(f"""
Dataset: Kuppe et al. (2022) Nature
- Total cells: {adata.n_obs:,}
- Total genes: {adata.n_vars:,}
- Samples: {adata.obs['sample'].nunique()}
- Donors: {adata.obs['donor_id'].nunique()}

Disease States:
- Myocardial infarction: {(adata.obs['disease'] == 'myocardial infarction').sum():,}
- Normal: {(adata.obs['disease'] == 'normal').sum():,}

Epicardial Markers Found: {len(markers)}/5
- {', '.join(markers.keys())}

Epicardial-related cells (by annotation): {epi_mask.sum():,}
- Note: Mostly adipocytes from epicardial fat
- True epicardial cells need marker-based identification
""")

    print("\nNext steps:")
    print("1. Use marker expression (WT1, TBX18, ALDH1A2) to identify epicardial cells")
    print("2. Focus on Fibroblast population which may contain epicardial cells")
    print("3. Integrate with Linna-Kuosmanen data (has mesothelial annotations)")

    print("\n" + "="*60)
    print("Analysis Complete!")
    print("="*60)

if __name__ == "__main__":
    main()
