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
| `communication_kuppe_sample.h5ad` | ✅ Created | Kuppe-only (1% sample) |
| `communication_merged_sample.h5ad` | ⚠️ Partial | Merged dataset, no batch correction yet |

### Full Data Files
| File | Status | Description |
|------|--------|-------------|
| `communication_kuppe.h5ad` | ✅ Created | Kuppe-only (175k cells) |
| `communication_merged.h5ad` | ⚠️ Created | Merged (757k cells), no batch correction |
| `epicardial_with_states.h5ad` | ✅ Created | 59k epicardial cells with state labels |

---

## Current State

### Phase 2: Data Preparation ⚠️ Almost Complete

**Completed:**
- Cell state classification (activated/quiescent/other)
- Communication dataset creation (kuppe + merged)
- Sample mode implemented (`--sample 0.01`)

**Issue: Harmony Batch Correction**

The Harmony integration is failing when storing results:
```
ValueError: Value passed for key 'X_pca_harmony' is of incorrect shape.
Values of obsm must match dimensions ('obs',) of parent.
Value had shape (50,) while it should have had (7571,).
```

**Debugging needed:**
```bash
python scripts/02_cell_states/prepare_communication_data.py --sample 0.01 --only merged --force
```

The script now prints debug info about `ho.Z_corr` shape. Need to check:
- `ho.Z_corr.shape` - expected (50, 7571) or (7571, 50)?
- May need `.cpu().numpy()` if it's a PyTorch tensor on MPS

**Workaround if needed:**
For development, batch correction is not critical. The sample data can be used without it.

---

### Phase 3: Communication Analysis 📝 Scripts Ready

Scripts created in `scripts/03_communication/`:

| Script | Purpose | Command |
|--------|---------|---------|
| `01_deg_analysis.py` | DEG with pyDESeq2 | `python 01_deg_analysis.py --sample` |
| `02_liana_analysis.py` | LIANA L-R inference | `python 02_liana_analysis.py --sample` |
| `03_nichenet_analysis.py` | NicheNet ligand activity | `python 03_nichenet_analysis.py --sample` |

**Dependencies needed:**
```bash
pip install pydeseq2 liana
```

**NicheNet data needed:**
Download from https://zenodo.org/record/7074291:
- `ligand_target_matrix.rds` → convert to CSV
- `lr_network.rds` → convert to CSV
Place in `data/nichenet/`

---

## Execution Order

### Step 1: Fix Harmony (or skip for now)
```bash
# Debug run
python scripts/02_cell_states/prepare_communication_data.py --sample 0.01 --only merged --force

# If Harmony keeps failing, the sample data is still usable for Phase 3 testing
```

### Step 2: DEG Analysis
```bash
pip install pydeseq2
python scripts/03_communication/01_deg_analysis.py --sample
# Output: results/deg/upregulated_genes_sample.txt
```

### Step 3: LIANA Analysis
```bash
pip install liana
python scripts/03_communication/02_liana_analysis.py --sample
# Output: results/liana/liana_epicardial_results_merged_sample.csv
```

### Step 4: NicheNet Analysis
```bash
# First, prepare NicheNet data (see above)
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

### 1. Numpy/scikit-misc Incompatibility
**Error:** `numpy.dtype size changed, may indicate binary incompatibility`
**Fix:** `pip install numpy==1.24.3 --force-reinstall`

### 2. HVG seurat_v3 Warning
**Warning:** `flavor='seurat_v3' expects raw count data, but non-integers were found`
**Status:** Safe to ignore - data is log-normalized, method still works

### 3. Harmony Output Shape
**Error:** `Value passed for key 'X_pca_harmony' is of incorrect shape`
**Status:** Debugging in progress - may need to handle PyTorch tensor conversion

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
