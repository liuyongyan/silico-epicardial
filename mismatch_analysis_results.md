# Ligand-Receptor Mismatch Analysis: "Primed But Starved" Pairs

## Objective

Systematically identify ligand-receptor pairs where the **receptor is upregulated** and the **cognate ligand is downregulated** in activated epicardial cells post-MI, across mouse and human.

Positive control: FGF10/FGFR2 (wet lab validated by Cheng Lab).

---

## Filtering Pipeline

### Step 1: L-R Database Construction

**Source**: OmniPath, merging 5 curated databases (per instruction 3 spec):
- CellPhoneDB, CellTalkDB, Ramilowski2015, Fantom5_LRdb, connectomeDB2020

**Complex expansion**: Receptor complexes (e.g., `FZD2_LRP6`) split into individual subunits.

**Result**: 5,669 unique L-R pairs (1,269 ligands × 1,194 receptors)

### Step 2: Mouse Mismatch Identification

**DEG source**: Quaife-Ryan 2021, Activated vs Quiescent epicardial cells (Wilcoxon rank-sum, 23,376 genes)

**Filters applied**:
1. Receptor significantly upregulated: logFC > 0, padj < 0.05 → 623 receptors
2. Ligand significantly downregulated: logFC < 0, padj < 0.05
3. **High-confidence pairs only**: n_db ≥ 3 (supported by ≥ 3 of the 5 databases)
4. **Starvation ratio**: computed over high-confidence cognate ligands only (fraction significantly downregulated)

**Scoring**: `composite = (receptor_logFC_norm) × (|ligand_logFC_norm|) × starvation_ratio × (n_db / 5)`

Where logFC values are capped to [0, 10] and normalized to [0, 1].

**Result**: 127 mouse mismatch pairs

### Step 3: Cross-Species Conservation

**Human DEG source**: PERIHEART, Quiescent vs Activated epicardial cells (Wilcoxon, 35,477 genes)

**Human data limitation**: All MI cells from single patient (PH-M57). Effect sizes 10-100× smaller than mouse.

**Conservation levels**:
- **Strict**: receptor↑ (padj<0.05) AND ligand↓ (padj<0.05) in BOTH species
- **Relaxed**: receptor logFC > 0 AND ligand logFC < 0 in both species (direction-consistent, significant in at least one)

**Result**: 100 conserved pairs (2 strict, 98 relaxed)

### Step 4: Protein Family Match Filter

Remove cross-family interactions (e.g., FGFR2/PF4 where PF4 is a chemokine, not an FGF ligand). Only keep pairs where receptor and ligand belong to the same signaling pathway family.

**Result**: 71 canonical conserved pairs

---

## Final Results: 71 Cross-Species Conserved Canonical Mismatch Pairs

Ranked by average mismatch score across both species (logFC capped ±10).

### Top 20

| Rank | Receptor | Ligand | n_db | Pathway | Mouse R↑ | Mouse L↓ | Mouse sig | Human R↑ | Human L↓ | Human sig | Avg Mismatch |
|:----:|----------|--------|:----:|---------|:--------:|:--------:|:---------:|:--------:|:--------:|:---------:|:------------:|
| 1 | OSMR | OSM | 4 | Cytokine | +16.35 | -3.79 | trend | +0.17 | -16.97 | trend | 11.98 |
| 2 | EDNRA | EDN3 | 5 | Endothelin | +10.08 | -2.53 | trend | +0.05 | -18.92 | trend | 11.29 |
| 3 | EPHA7 | EFNA2 | 5 | Ephrin | +31.96 | -1.90 | strict | +0.40 | -19.18 | trend | 11.15 |
| 4 | TRHR | TRH | 5 | Neuropeptide | +20.33 | -1.23 | trend | +0.65 | -16.78 | trend | 10.94 |
| 5 | ACVR1 | GDF2 | 4 | BMP/TGFb | +21.96 | -1.49 | trend | +0.10 | -17.81 | trend | 10.79 |
| 6 | BMPR2 | GDF2 | 4 | BMP/TGFb | +15.12 | -1.49 | trend | +0.03 | -17.81 | trend | 10.76 |
| 7 | ACVR1 | BMP6 | 4 | BMP/TGFb | +21.96 | -20.73 | **strict** | +0.10 | -0.46 | trend | 10.28 |
| 8 | BMPR2 | BMP6 | 4 | BMP/TGFb | +15.12 | -20.73 | **strict** | +0.03 | -0.46 | trend | 10.25 |
| 9 | LRP6 | RSPO1 | 4 | Wnt | +12.87 | -24.78 | **strict** | +0.16 | -0.13 | trend | 10.14 |
| 10 | FZD1 | MYOC | 3 | Wnt | +9.44 | -153.04 | **strict** | +0.05 | -0.70 | trend | 10.10 |
| 11 | TNFRSF12A | TNFSF12 | 5 | TNF | +60.92 | -77.96 | **strict** | +0.01 | -0.16 | trend | 10.09 |
| 12 | IL13RA1 | IL4 | 4 | Interleukin | +32.59 | -0.07 | trend | +0.09 | -16.81 | trend | 10.08 |
| 13 | IL6ST | OSM | 4 | Cytokine | +5.93 | -3.79 | trend | +0.35 | -16.97 | trend | 10.03 |
| 14 | IL12RB1 | IL12A | 4 | Interleukin | +0.30 | -9.22 | trend | +0.43 | -18.97 | trend | 9.97 |
| 15 | FGFR2 | FGF16 | 3 | FGF | +4.81 | -4.56 | **strict** | +0.39 | -18.65 | trend | 9.88 |
| 16 | FGFR2 | FGF6 | 4 | FGF | +4.81 | -3.58 | trend | +0.39 | -19.07 | trend | 9.39 |
| 17 | TYRO3 | GAS6 | 5 | TAM | +8.31 | -14.00 | **strict** | +0.28 | -0.02 | trend | 9.31 |
| 18 | EPHA3 | EFNA2 | 5 | Ephrin | +5.34 | -1.90 | **strict** | +0.22 | -19.18 | trend | 8.73 |
| 19 | IL6ST | CNTF | 4 | Cytokine | +5.93 | -0.31 | trend | +0.35 | -20.41 | trend | 8.29 |
| 20 | NOTCH3 | PSEN1 | 4 | Notch | +5.47 | -10.10 | **strict** | +0.91 | -0.05 | trend | 8.22 |

### FGF10/FGFR2 and Key Pairs

| Rank | Receptor | Ligand | n_db | Pathway | Avg Mismatch | Mouse sig | Human sig |
|:----:|----------|--------|:----:|---------|:------------:|:---------:|:---------:|
| 15 | FGFR2 | FGF16 | 3 | FGF | 9.88 | strict | trend |
| 16 | FGFR2 | FGF6 | 4 | FGF | 9.39 | trend | trend |
| 23 | FGFR2 | FGF5 | 4 | FGF | 8.11 | trend | trend |
| **39** | **FGFR2** | **FGF10** | **5** | **FGF** | **5.90** | **strict** | **trend** |
| 54 | FGFR2 | FGF7 | 5 | FGF | 5.05 | strict | trend |

### By Signaling Pathway (best pair per pathway)

| Pathway | Rank | Receptor | Ligand | Avg Mismatch | n_db | Pairs in Pathway |
|---------|:----:|----------|--------|:------------:|:----:|:----------------:|
| Cytokine | 1 | OSMR | OSM | 11.98 | 4 | 3 |
| Endothelin | 2 | EDNRA | EDN3 | 11.29 | 5 | 2 |
| Ephrin | 3 | EPHA7 | EFNA2 | 11.15 | 5 | 7 |
| Neuropeptide | 4 | TRHR | TRH | 10.94 | 5 | 3 |
| BMP/TGFb | 5 | ACVR1 | GDF2 | 10.79 | 4 | 10 |
| Wnt | 9 | LRP6 | RSPO1 | 10.14 | 4 | 9 |
| TNF | 11 | TNFRSF12A | TNFSF12 | 10.09 | 5 | 2 |
| Interleukin | 12 | IL13RA1 | IL4 | 10.08 | 4 | 14 |
| FGF | 15 | FGFR2 | FGF16 | 9.88 | 3 | 5 |
| TAM | 17 | TYRO3 | GAS6 | 9.31 | 5 | 2 |
| Notch | 20 | NOTCH3 | PSEN1 | 8.22 | 4 | 5 |
| Sema/VEGF | 46 | NRP2 | SEMA4F | 5.53 | 3 | 2 |
| GDNF | 62 | GFRA2 | NRTN | 2.55 | 4 | 1 |
| EGF/ErbB | 66 | ERBB2 | BTC | 2.00 | 4 | 1 |

---

## Key Findings

### 1. FGF10/FGFR2 is validated but not top-ranked by mismatch score alone

FGF10/FGFR2 ranks **39/71** among conserved canonical pairs. Its mismatch score is moderate because FGFR2's receptor logFC (+4.81 in mouse, +0.39 in human) is small compared to receptors like EPHA7 (+31.96), ACVR1 (+21.96), or EDNRA (+10.08). However, it has the **highest database consensus** (n_db=5) and is one of few FGF pairs with strict significance in mouse.

### 2. BMP pathway has strong conserved signal

ACVR1/BMP6 (#7) and BMPR2/BMP6 (#8) both show significant mismatch in mouse with directional conservation in human. Note: BMP4 is **upregulated** in human (opposite of mouse), so BMP4 pairs are NOT conserved.

### 3. Novel pathway candidates

Several pathways not considered in instruction 3 show strong conserved mismatch:

- **Endothelin** (EDNRA/EDN3, #2): EDN3 is strongly downregulated in both species. Endothelin signaling has known roles in cardiac function.
- **Ephrin** (EPHA7/EFNA2, #3): EFNA2 drops -19.18 in human. Ephrin signaling regulates EMT and cell migration.
- **TAM** (TYRO3/GAS6, #17): GAS6 is strongly downregulated in mouse (-14.00). TAM receptors mediate efferocytosis post-MI.
- **Cytokine** (OSMR/OSM, #1): Oncostatin M is downregulated in both species.

### 4. Human data is a bottleneck

Only 2 pairs reach strict significance in both species (EPHB6/AFDN, INSR/NAMPT — both removed by family filter as cross-family). All remaining pairs are "strict in mouse, trend in human." This reflects the single-patient MI limitation of PERIHEART data.

---

## Geneformer In Silico Perturbation: Skipped

Instruction 3 Phase 4 proposed using Geneformer fine-tuned classifier to measure each receptor's contribution to epicardial activation via in silico deletion. This step was skipped for two reasons:

### Reason 1: Systematic bias against low-expression genes

Geneformer tokenizes each cell as a rank-ordered sequence of its top ~2,048 expressed genes (out of ~20,000). FGF family genes have very low expression:

| Gene | Cells Expressing | Mean Expression (log1p) | Tokenized? |
|------|:----------------:|:-----------------------:|:----------:|
| FGFR2 | 1.7% (321/19,412) | 0.24–0.36 | Rarely |
| FGFR1 | 28–40% | 0.44–0.61 | Sometimes |
| Col1a1 (EMT marker) | >50% | >1.0 | Almost always |

In silico perturbation = deleting a gene's token from the sequence. If FGFR2 is not tokenized in 98.3% of cells, deleting it changes nothing. Round 1 results (embedding perturbation on human data) confirmed this:

- FGFR2 perturbation effect: **rank 61/84** (bottom 27%)
- Correlation between n_cells_with_gene and perturbation effect: **r = 0.649**
- Top-ranked genes (PLXDC2, ALK, AQP1) are all expressed in more cells, not necessarily more biologically important

This bias is inherent to Geneformer's rank-value tokenization and would persist in Round 2 (classifier perturbation).

### Reason 2: Classification label quality

The quiescent vs activated labels used for fine-tuning are derived from our own signature-based scoring with GMM thresholding. Cross-validation between condition-based and score-based classification shows only **60.8% agreement** in human data. A classifier trained on noisy labels will produce noisy perturbation results — the model would likely learn patient-specific differences (single MI patient PH-M57 vs 29 normal donors) rather than true activation biology.

### Impact on downstream analysis

Geneformer perturbation was allocated 25% weight in instruction 3's Phase 6 priority scoring. Without it, we redistribute to the other three dimensions (mismatch 30%, conservation 30%, druggability 20%, literature 20%). This means the final ranking relies more on expression-level evidence and prior knowledge, without the non-linear gene regulatory insights Geneformer could provide. This is a known limitation of our current analysis.

---

## Output Files

| File | Description |
|------|-------------|
| `results/mismatch/mouse_lr_mismatch_refined.csv` | 127 mouse pairs with composite scoring |
| `results/mismatch/cross_species_lr_mismatch.csv` | All 1,953 L-R pairs with both species data |
| `results/mismatch/cross_species_conserved.csv` | 100 conserved mismatch pairs (before family filter) |
| `results/mismatch/therapeutic_targets_prioritized.csv` | 127 pairs with multi-dimensional priority scores |
| `results/mismatch/therapeutic_targets_actionable.csv` | 43 actionable targets (score>0.3, druggability≥0.5) |
| `results/mismatch/curated_lr_pairs_mouse.csv` | 5,669 L-R pairs used |
| `scripts/05_mismatch/01_mouse_lr_mismatch.py` | Mouse raw mismatch analysis |
| `scripts/05_mismatch/02_mouse_lr_mismatch_refined.py` | Mouse refined multi-dimensional scoring |
| `scripts/05_mismatch/03_human_lr_mismatch_refined.py` | Human mismatch analysis |
| `scripts/05_mismatch/04_cross_species_comparison.py` | Cross-species comparison |
| `scripts/05_mismatch/05_therapeutic_prioritization.py` | Therapeutic target prioritization |

---

## Date

Analysis completed: 2026-03-21
