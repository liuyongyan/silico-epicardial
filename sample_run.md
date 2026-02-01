# Sample Run Strategy & Current State

## Overview

This document tracks the development workflow using **1% sampled data** to quickly iterate through the pipeline before running with full data.

## Sample Strategy

### Why Sample?
- Full dataset: ~757k cells → 30-80 min per script
- 1% sample: ~7.5k cells → 3-5 min per script
- Goal: Get all scripts working end-to-end, then run full data once

### Sample Files
| File | Status | Description |
|------|--------|-------------|
| `communication_kuppe_sample.h5ad` | ⚠️ Needs regen | Kuppe-only (1% sample) - needs disease filtering |
| `communication_merged_sample.h5ad` | ⚠️ Needs regen | Merged dataset - needs disease filtering |

### Full Data Files
| File | Status | Description |
|------|--------|-------------|
| `epicardial_periheart.h5ad` | ⚠️ Needs regen | PERIHEART epicardial - needs disease filtering |
| `epicardial_carebank.h5ad` | ⚠️ Needs regen | CAREBANK epicardial - needs disease filtering |
| `epicardial_kuppe.h5ad` | ⚠️ Needs regen | Kuppe epicardial - needs disease filtering |
| `epicardial_with_states.h5ad` | ⚠️ Needs regen | All epicardial cells with state labels |
| `communication_kuppe.h5ad` | ⚠️ Needs regen | Kuppe-only communication dataset |
| `communication_merged.h5ad` | ⚠️ Needs regen | Merged communication dataset |

---

## Current State

### Phase 1: Preprocessing ⚠️ Needs Re-run

**Status:** Disease filtering added - need to re-run extraction scripts.

Scripts modified:
- `extract_epicardial_simple.py` - now filters to MI + Normal
- `extract_kuppe_epicardial.py` - now filters to MI + Normal

---

### Phase 2: Data Preparation ⚠️ Needs Re-run

**Completed:**
- Cell state classification (activated/quiescent/other)
- Communication dataset creation (kuppe + merged)
- Sample mode implemented (`--sample 0.01`)
- Harmony batch correction working
- **Disease filtering added to sender cell extraction**

**Harmony Fix (Jan 2026):**
The original code incorrectly transposed `ho.Z_corr`. harmonypy returns `(n_cells, n_pcs)` directly, not `(n_pcs, n_cells)` as the comment suggested. Removed the `.T` transpose to fix.

**Disease Filtering (Jan 2026):**
Added disease filtering to `prepare_communication_data.py` to keep only MI + Normal sender cells.

---

### Phase 3: Communication Analysis ✅ Mostly Complete

Scripts in `scripts/03_communication/`:

| Script | Status | Key Results |
|--------|--------|-------------|
| `01_deg_analysis.py` | ✅ Done | 439 upregulated, 2,714 downregulated genes |
| `02_liana_analysis.py` | ✅ Done | 2,104 L-R pairs targeting epicardial cells |
| `03_nichenet_analysis.py` | ⏸️ Needs data | Requires NicheNet data from Zenodo |

**DEG Top Hits (Upregulated in Activated):**
- EMT/ECM: LAMA2, COL6A3, COL4A1/A2, FN1, VCAN, FBN1
- EMT TFs: ZEB1, ZEB2 significantly upregulated

**LIANA Top L-R Pairs (→ Epicardial):**
- FGF12 → FGFR1/FGFR2 (Cardiomyocyte)
- TGFB1 → ACVR1B_TGFBR2 (Epicardial autocrine)
- NRG1/NRG3 → EGFR (Endothelial)
- DCN → EGFR/ERBB4/MET (Fibroblast)
- Top ligands: CALM1, FN1, TGFB1, LAMA2, ADAM10

**NicheNet data (if needed):**
Download from https://zenodo.org/record/7074291:
- `ligand_target_matrix.rds` → convert to CSV
- `lr_network.rds` → convert to CSV
Place in `data/nichenet/`

---

## Execution Order

### Step 0: Re-run Preprocessing with Disease Filtering
```bash
# Disease filtering added to keep only MI + Normal
# Re-extract epicardial cells (Phase 1)
python scripts/01_preprocessing/extract_epicardial_simple.py   # PERIHEART + CAREBANK
python scripts/01_preprocessing/extract_kuppe_epicardial.py    # Kuppe

# Re-classify cell states (Phase 2)
python scripts/02_cell_states/classify_cell_states.py
```

### Step 1: ✅ Harmony Fixed
```bash
# Sample data now includes batch correction + disease filtering
python scripts/02_cell_states/prepare_communication_data.py --sample 0.01 --only merged --force
```

### Step 2: ✅ DEG Analysis Complete
```bash
python scripts/03_communication/01_deg_analysis.py --sample
# Output: results/deg/deg_results_sample.csv (gene symbols)
```

### Step 3: ✅ LIANA Analysis Complete
```bash
python scripts/03_communication/02_liana_analysis.py --sample
# Output: results/liana/liana_epicardial_results_merged_sample.csv
```

### Step 4: NicheNet Analysis (Optional)
```bash
# First, prepare NicheNet data from Zenodo
python scripts/03_communication/03_nichenet_analysis.py --sample
# Output: results/nichenet/nichenet_ligand_activities_sample.csv
```

### Step 5: Run Full Pipeline
After all scripts work with sample data:
```bash
# Regenerate full communication_merged.h5ad with Harmony
python scripts/02_cell_states/prepare_communication_data.py --force

# Run Phase 3 without --sample flag
python scripts/03_communication/01_deg_analysis.py
python scripts/03_communication/02_liana_analysis.py
python scripts/03_communication/03_nichenet_analysis.py
```

---

## Known Issues & Fixes

### 0. Disease Condition Filtering (Jan 2026)
**Issue:** CAREBANK dataset contains no MI patients (only myocardial ischemia, heart valve disorder, etc.), confounding DEG and LIANA results.
**Fix:** Added disease filtering to all extraction scripts:
- `extract_epicardial_simple.py` - filters PERIHEART/CAREBANK
- `extract_kuppe_epicardial.py` - filters Kuppe
- `prepare_communication_data.py` - filters sender cells
Only keeps `['myocardial infarction', 'normal']` conditions.

### 1. Numpy/scikit-misc Incompatibility
**Error:** `numpy.dtype size changed, may indicate binary incompatibility`
**Fix:** `pip install numpy==1.24.3 --force-reinstall`

### 2. HVG seurat_v3 Warning
**Warning:** `flavor='seurat_v3' expects raw count data, but non-integers were found`
**Status:** Safe to ignore - data is log-normalized, method still works

### 3. Harmony Output Shape
**Error:** `Value passed for key 'X_pca_harmony' is of incorrect shape`
**Status:** ✅ Fixed - removed incorrect `.T` transpose; `ho.Z_corr` is already `(n_cells, n_pcs)`

---

## File Structure

```
Silico Epicardial/
├── data/
│   ├── processed/
│   │   ├── epicardial_with_states.h5ad      # Phase 2 output
│   │   ├── communication_kuppe.h5ad         # Full Kuppe
│   │   ├── communication_kuppe_sample.h5ad  # 1% Kuppe
│   │   ├── communication_merged.h5ad        # Full merged
│   │   └── communication_merged_sample.h5ad # 1% merged
│   └── nichenet/                            # NicheNet data (to download)
├── scripts/
│   ├── 02_cell_states/
│   │   ├── classify_cell_states.py
│   │   └── prepare_communication_data.py
│   └── 03_communication/
│       ├── 01_deg_analysis.py
│       ├── 02_liana_analysis.py
│       └── 03_nichenet_analysis.py
└── results/
    ├── deg/
    ├── liana/
    └── nichenet/
```

---

## Quick Reference

```bash
# Check sample data
python -c "import anndata as ad; a=ad.read_h5ad('data/processed/communication_merged_sample.h5ad'); print(a)"

# Run all Phase 3 with sample
for script in scripts/03_communication/*.py; do python "$script" --sample; done
```
