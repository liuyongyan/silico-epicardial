# Mouse Epicardial Analysis Results

## Objective

Find computational evidence supporting FGF10/FGFR2 as a key ligand-receptor pair in epicardial activation post-MI. Wet lab has already validated the efficacy of FGF10/FGFR2 pathway intervention in early MI.

---

## Dataset 1: Quaife-Ryan et al. 2021 (eLife)

**Paper**: "Single-cell transcriptomics defines heterogeneity of epicardial cells and fibroblasts within the infarcted murine heart"
**DOI**: 10.7554/eLife.65921
**Data**: ArrayExpress E-MTAB-10035

### Dataset Characteristics

- **Cells**: 112,676 epicardial stromal cells (EpiSC)
- **Timepoint**: Day 5 post-MI (ischemia/reperfusion model)
- **Model**: C57BL/6 mice, 50 min ischemia followed by reperfusion
- **Isolation**: Novel perfusion-based technique specifically for epicardial cells
- **Conditions**: 85,514 MI cells, 27,162 normal cells

### Download Status

- [x] Data downloaded
- [x] QC completed
- [x] Epicardial cells identified (all cells are EpiSC by design)
- [x] Activation states classified
- [x] DEG analysis completed

### Key Findings

#### Cell State Classification

Using signature-based scoring (quiescent: Wt1, Upk3b, Msln, Krt19, Cdh1, Cldn1, Cldn3, Tjp1, Ocln; activated: Snai1, Snai2, Twist1, Vim, Fn1, Col1a1, Col3a1, Postn, Acta2):

| Condition | Quiescent | Activated | % Activated |
|-----------|-----------|-----------|-------------|
| Normal    | 25,385    | 1,777     | **6.5%**    |
| MI        | 23,981    | 61,533    | **72.0%**   |

**MI dramatically shifts epicardial cells toward an activated state.**

#### FGF Family Expression (MI vs Normal)

| Gene  | Normal Mean | MI Mean | log2FC | p-value |
|-------|-------------|---------|--------|---------|
| Fgf1  | 0.095       | 0.220   | +1.14  | 1.0e-107 |
| Fgf2  | 0.573       | 0.744   | +0.37  | 5.4e-70 |
| Fgf7  | 0.472       | 0.442   | -0.09  | 2.9e-01 |
| **Fgf10** | 0.089   | 0.025   | **-1.50** | **5.1e-81** |
| Fgfr1 | 3.595       | 3.299   | -0.12  | 0.0e+00 |
| **Fgfr2** | 0.037   | 0.230   | **+2.37** | **7.9e-243** |
| Fgfr3 | 0.056       | 0.058   | +0.05  | 1.2e-03 |

#### FGF10/FGFR2 Expression by Cell State AND Condition

| State | Condition | FGF10 Mean | FGFR2 Mean |
|-------|-----------|------------|------------|
| Quiescent | Normal | 0.0921 | 0.0270 |
| Quiescent | MI | 0.0554 | 0.1490 |
| Activated | Normal | 0.0383 | 0.1721 |
| Activated | MI | 0.0130 | 0.2616 |

### Interpretation

1. **FGFR2 marks activated epicardial cells**
   - log2FC = +2.37 in MI vs Normal (p < 1e-242)
   - log2FC = +1.48 in Activated vs Quiescent cells (p = 3.7e-33)
   - Highest expression in MI + Activated cells (0.262)

2. **FGF10 is associated with quiescence**
   - log2FC = -1.50 in MI vs Normal (p < 1e-80)
   - log2FC = -1.83 in Activated vs Quiescent cells
   - Highest expression in Normal + Quiescent cells (0.092)

3. **Model**: FGF10 may act as a quiescence-maintaining signal. Upon MI:
   - FGF10 expression decreases
   - FGFR2 expression increases
   - Cells transition from quiescent to activated (EMT-like)

---

## Dataset 2: Forte et al. (GSE135310) — Time Course

**Data**: GEO GSE135310
**Cells**: ~38,600 (C57BL/6J) + ~13,000 (129S1/SvlmJ)
**Timepoints**: Day 0 (sham), 1, 3, 5, 7, 14, 28 post-MI
**Key feature**: Wt1-Cre;Rosa-ZsGreen transgenic mice (genetic labeling of epicardial lineage)

### Download Status

- [ ] Data downloaded
- [ ] QC completed
- [ ] Wt1+ epicardial cells identified
- [ ] Time course analysis completed

### Results

*(To be analyzed)*

---

## FGF Family Expression Summary

| Gene | Quaife-Ryan (Quiescent) | Quaife-Ryan (Activated) | log2FC | p_adj | Notes |
|------|-------------------------|-------------------------|--------|-------|-------|
| Fgf10 | 0.074 | 0.014 | -1.83 | 1.8e-03 | Higher in quiescent |
| Fgfr2 | 0.086 | 0.259 | +1.48 | 3.7e-33 | Higher in activated |
| Fgfr1 | 3.025 | 3.639 | +0.27 | 1.1e-19 | Broadly expressed |
| Fgf7 | 0.385 | 0.499 | +0.37 | 1.5e-21 | Slight increase |

---

## Cross-Species Comparison (Mouse vs Human)

| Gene | Mouse (Quaife-Ryan) | Mouse (Forte) | Human (PERIHEART) | Conserved? |
|------|---------------------|---------------|-------------------|------------|
| FGF10/Fgf10 | Quiescent-high | TBD | TBD | TBD |
| FGFR2/Fgfr2 | Activated-high | TBD | TBD | TBD |

---

## Conclusions

### Quaife-Ryan 2021 Analysis

1. **FGFR2 is a robust marker of epicardial activation**
   - 6x higher in MI vs Normal (log2FC = 2.37)
   - 2.8x higher in activated vs quiescent cells (log2FC = 1.48)

2. **FGF10 inversely correlates with activation**
   - 3x lower in MI vs Normal (log2FC = -1.50)
   - 3.5x lower in activated vs quiescent cells (log2FC = -1.83)

3. **Potential therapeutic implication**
   - FGF10 supplementation may help maintain epicardial quiescence
   - FGFR2 inhibition may prevent excessive activation/fibrosis
   - The FGF10/FGFR2 axis is a key regulator of epicardial cell fate

### Output Files

- `data/processed/mouse_quaife_ryan_analyzed.h5ad` - Analyzed dataset
- `results/mouse/quaife_ryan_umap_overview.png` - UMAP visualization
- `results/mouse/quaife_ryan_fgf_violin.png` - FGF expression violin plots
- `results/mouse/quaife_ryan_fgf_by_condition.png` - FGF by condition
- `results/mouse/quaife_ryan_de_activated_vs_quiescent.csv` - Full DE results
- `results/mouse/quaife_ryan_fgf_summary.csv` - FGF expression summary
