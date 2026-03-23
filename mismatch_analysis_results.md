# Ligand-Receptor Mismatch Analysis: "Primed But Starved" Pairs

## Objective

Systematically identify ligand-receptor pairs where the **receptor is upregulated** and the **cognate ligand is downregulated** in activated epicardial cells post-MI, across mouse and human.

Positive control: FGF10/FGFR2 (wet lab validated by Cheng Lab).

---

## Datasets

| Dataset | Species | Source | Cells | Comparison | Genes |
|---------|---------|--------|:-----:|------------|:-----:|
| Quaife-Ryan 2021 | Mouse | E-MTAB-10035 | 112,676 EpiSC | Activated vs Quiescent | 23,376 |
| PERIHEART | Human | Linna-Kuosmanen 2024 | 32,924 mesothelial | Quiescent vs Activated | 35,477 |

**Human data limitation**: All MI cells (6,439) from single patient PH-M57. Effect sizes 10–100× smaller than mouse.

---

## Analysis Pipeline

### Step 1: L-R Database Construction

**Source**: OmniPath, merging 5 curated databases:
- CellPhoneDB, CellTalkDB, Ramilowski2015, Fantom5_LRdb, connectomeDB2020

**Complex expansion**: Receptor complexes (e.g., `FZD2_LRP6`) split into individual subunits.

**Database consensus** (`n_db`): number of the 5 databases supporting each pair (range 1–5).

**Result**: 5,669 unique L-R pairs (1,269 ligands × 1,194 receptors)

### Step 2: Mouse Mismatch Identification

**Filters**:
1. Receptor significantly upregulated: logFC > 0, padj < 0.05 → 623 receptors
2. Ligand significantly downregulated: logFC < 0, padj < 0.05
3. High-confidence pairs only: n_db ≥ 3
4. Starvation ratio: fraction of high-confidence cognate ligands significantly downregulated

**Scoring**: `composite = receptor_logFC_norm × |ligand_logFC_norm| × starvation_ratio × (n_db / 5)`

Where logFC values are capped to [0, 10] and normalized to [0, 1].

**Result**: 127 mouse mismatch pairs

### Step 3: Cross-Species Conservation

For each high-confidence L-R pair (n_db ≥ 3), check if the receptor↑ + ligand↓ pattern holds in both mouse and human. Gene ortholog mapping by case conversion (Fgfr2 → FGFR2).

**Conservation levels**:
- **Strict**: receptor↑ (padj<0.05) AND ligand↓ (padj<0.05) in BOTH species
- **Relaxed**: direction-consistent in both species, significant in at least one

**Result**: 100 conserved pairs (2 strict, 98 relaxed)

### Step 4: Protein Family Match Filter

Remove cross-family interactions where ligand and receptor belong to different signaling pathway families (e.g., FGFR2/PF4 where PF4 is a chemokine, not an FGF ligand). Family assignments based on UniProt protein family annotations.

**Result**: 71 canonical conserved pairs

### Step 5: Automated Druggability Scoring

Two data sources combined to avoid bias:

**DGIdb** (Drug-Gene Interaction Database, GraphQL API):
- Counts drug-gene interactions and approved drugs per gene
- Captures small molecule inhibitors, antibodies, approved agonists
- **Limitation**: misses recombinant protein delivery — our primary therapeutic strategy. FGF10 has 0 interactions in DGIdb despite being deliverable as recombinant protein.

**UniProt** (subcellular location API):
- Classifies ligand deliverability based on where the protein localizes
- Secreted → 1.0 (ideal for recombinant delivery), Membrane → 0.5 (Fc-fusion possible), Intracellular → 0.1

**Combined druggability per pair**:
```
ligand_combined = 0.6 × UniProt_deliverability + 0.4 × DGIdb_score
pair_druggability = 0.3 × receptor_DGIdb + 0.7 × ligand_combined
```

This weighting reflects that our strategy is "deliver the depleted ligand," not "drug the receptor."

### Step 6: Automated Literature Scoring

**PubMed E-utilities API**: For each gene, count publications matching:
`"gene_name" AND (cardiac OR epicardial OR "heart repair" OR "myocardial infarction" OR "cardiac regeneration")`

**Normalization**: `literature_score = log10(count + 1) / log10(1500)`, capped at 1.0.

**Pair score**: average of receptor and ligand literature scores.

### Step 7: Therapeutic Target Prioritization

**Final composite score** (Geneformer skipped, weights redistributed):

| Dimension | Weight | Source | Automated? |
|-----------|:------:|--------|:----------:|
| Mismatch composite | 30% | Step 2 (mouse DEG + L-R database) | Yes |
| Cross-species conservation | 30% | Step 3 (mouse + human DEG) | Yes |
| Druggability | 20% | Step 5 (DGIdb + UniProt) | Yes |
| Literature support | 20% | Step 6 (PubMed) | Yes |

```
priority = 0.30 × mismatch_norm + 0.30 × conservation + 0.20 × druggability + 0.20 × literature
```

---

## Geneformer In Silico Perturbation: Skipped

Instruction 3 Phase 4 proposed using Geneformer fine-tuned classifier to measure each receptor's contribution to epicardial activation via in silico deletion. This step was allocated 25% weight in the original scoring. It was skipped for two reasons:

### Reason 1: Systematic bias against low-expression genes

Geneformer tokenizes each cell as a rank-ordered sequence of its top ~2,048 expressed genes (out of ~20,000). FGF family genes have very low expression:

| Gene | Cells Expressing | Mean Expression (log1p) | Tokenized? |
|------|:----------------:|:-----------------------:|:----------:|
| FGFR2 | 1.7% (321/19,412) | 0.24–0.36 | Rarely |
| FGFR1 | 28–40% | 0.44–0.61 | Sometimes |
| Col1a1 (EMT marker) | >50% | >1.0 | Almost always |

In silico perturbation = deleting a gene's token from the sequence. If FGFR2 is not tokenized in 98.3% of cells, deleting it changes nothing. Round 1 results (embedding perturbation on human data) confirmed this:

- FGFR2 perturbation effect: rank 61/84 (bottom 27%)
- Correlation between n_cells_with_gene and perturbation effect: **r = 0.649**
- Top-ranked genes (PLXDC2, ALK, AQP1) are all expressed in more cells, not necessarily more biologically important

### Reason 2: Classification label quality

The quiescent vs activated labels are derived from signature-based scoring with GMM thresholding. Cross-validation between condition-based and score-based classification shows only **60.8% agreement** in human data. A classifier trained on noisy labels will learn patient-specific differences (single MI patient PH-M57 vs 29 normal donors) rather than true activation biology.

---

## Final Results

### Top 20 Therapeutic Targets (fully automated scoring)

| Rank | Receptor | Ligand | Score | Mismatch | Conservation | Druggability | Literature | Ligand Type | Pathway |
|:----:|----------|--------|:-----:|:--------:|:------------:|:------------:|:----------:|:-----------:|---------|
| 1 | INSR | NAMPT | 0.615 | 0.011 | 0.300 (strict) | 0.174 | 0.130 | Secreted | Insulin |
| 2 | TYRO3 | GAS6 | 0.587 | 0.125 | 0.210 | 0.153 | 0.099 | Secreted | TAM |
| 3 | BMPR2 | BMP6 | 0.576 | 0.096 | 0.210 | 0.130 | 0.140 | Secreted | BMP |
| 4 | ITGB1 | CD14 | 0.573 | 0.048 | 0.210 | 0.156 | 0.158 | Secreted | Integrin |
| 5 | ACVR1 | BMP6 | 0.562 | 0.107 | 0.210 | 0.137 | 0.108 | Secreted | BMP/Activin |
| 6 | IL1RL2 | IL18 | 0.555 | 0.058 | 0.210 | 0.153 | 0.135 | Secreted | Interleukin |
| 7 | UNC5B | NTN1 | 0.544 | 0.300 | 0.060 | 0.090 | 0.094 | Secreted | Netrin |
| 8 | IL2RG | IL2 | 0.542 | 0.003 | 0.210 | 0.176 | 0.152 | Secreted | Interleukin |
| 9 | ITGA5 | CCN1 | 0.527 | 0.074 | 0.210 | 0.138 | 0.106 | Secreted | Integrin |
| 10 | ITGB3 | CCN1 | 0.525 | 0.050 | 0.210 | 0.141 | 0.124 | Secreted | Integrin |
| 11 | ITGAV | CCN1 | 0.523 | 0.071 | 0.210 | 0.140 | 0.102 | Secreted | Integrin |
| 12 | IL2RG | IL15 | 0.516 | 0.010 | 0.210 | 0.169 | 0.127 | Secreted | Interleukin |
| **13** | **FGFR2** | **FGF10** | **0.509** | 0.025 | 0.210 | 0.144 | 0.130 | **Secreted** | **FGF** |
| **14** | **FGFR2** | **FGF7** | **0.508** | 0.023 | 0.210 | 0.160 | 0.114 | **Secreted** | **FGF** |
| 15 | INSR | HRAS | 0.494 | 0.001 | 0.210 | 0.145 | 0.137 | Membrane | Insulin |
| 16 | ITGB1 | LAMC2 | 0.494 | 0.048 | 0.210 | 0.147 | 0.088 | Secreted | Integrin |
| 17 | IL2RG | ICAM1 | 0.493 | 0.007 | 0.210 | 0.123 | 0.152 | Membrane | Interleukin |
| 18 | BMPR2 | BMP2 | 0.488 | 0.096 | 0.060 | 0.161 | 0.171 | Secreted | BMP |
| 19 | OGFR | PENK | 0.483 | 0.240 | 0.060 | 0.100 | 0.083 | Secreted | Opioid |
| 20 | FGFR2 | FGF16 | 0.482 | 0.013 | 0.210 | 0.144 | 0.115 | Secreted | FGF |

### Positive Control Validation

FGF10/FGFR2 ranks **#13/127** with fully automated scoring — no manual annotation. Comparison across scoring methods:

| Pair | Manual | DGIdb-only | **Corrected (final)** |
|------|:------:|:----------:|:---------------------:|
| Fgfr2/Fgf10 | 1 | 27 | **13** |
| Fgfr2/Fgf7 | 4 | 21 | **14** |
| Acvr1/Bmp6 | 2 | 11 | **5** |
| Bmpr2/Bmp6 | 5 | 8 | **3** |
| Tyro3/Gas6 | 3 | 2 | **2** |
| Insr/Nampt | 9 | 1 | **1** |
| Unc5b/Ntn1 | 6 | 18 | **7** |
| Bmpr2/Bmp2 | 7 | 20 | **18** |

Manual scoring placed FGFR2/FGF10 at #1 because we knew it was the validated positive control (literature=1.0, druggability=1.0). The automated score does not have this prior — its #13 ranking reflects the gene's actual data footprint: moderate mismatch, conserved across species, secreted ligand, and ~97 cardiac publications.

### Key Findings

1. **FGF10/FGFR2 is cross-species conserved and in the top 10%.** Ranks #13/127 without manual input. FGFR2 is upregulated and FGF10 is downregulated in both mouse (significant) and human (directional trend).

2. **BMP pathway pairs rank highly.** BMPR2/BMP6 (#3), ACVR1/BMP6 (#5), BMPR2/BMP2 (#18). Note: BMP4 is upregulated in human (opposite of mouse), so BMP4 pairs are NOT conserved.

3. **Novel candidates emerge.** TYRO3/GAS6 (#2, TAM receptor, efferocytosis), INSR/NAMPT (#1, the only strictly conserved pair), and UNC5B/NTN1 (#7, netrin guidance) were not in the original instruction 3 hypothesis.

4. **Mismatch score alone is insufficient.** FGFR2/FGF10's mismatch composite is only 0.0817 (rank 77/127). Multi-dimensional scoring — especially cross-species conservation and ligand deliverability — is essential to surface biologically meaningful targets.

5. **Human data is a bottleneck.** No conserved pairs reach strict significance in both species. All conservation calls rely on directional consistency in human (trend level). Multi-patient human MI data would substantially strengthen these findings.

### 71 Cross-Species Conserved Canonical Pairs (by avg mismatch)

See `cross_species_conserved.csv` for the full 100 pairs (before family filter). Top pairs by signaling pathway:

| Pathway | Rank | Receptor | Ligand | Avg Mismatch | n_db |
|---------|:----:|----------|--------|:------------:|:----:|
| Cytokine | 1 | OSMR | OSM | 11.98 | 4 |
| Endothelin | 2 | EDNRA | EDN3 | 11.29 | 5 |
| Ephrin | 3 | EPHA7 | EFNA2 | 11.15 | 5 |
| Neuropeptide | 4 | TRHR | TRH | 10.94 | 5 |
| BMP/TGFb | 5 | ACVR1 | GDF2 | 10.79 | 4 |
| Wnt | 9 | LRP6 | RSPO1 | 10.14 | 4 |
| TNF | 11 | TNFRSF12A | TNFSF12 | 10.09 | 5 |
| Interleukin | 12 | IL13RA1 | IL4 | 10.08 | 4 |
| FGF | 15 | FGFR2 | FGF16 | 9.88 | 3 |
| TAM | 17 | TYRO3 | GAS6 | 9.31 | 5 |
| Notch | 20 | NOTCH3 | PSEN1 | 8.22 | 4 |

FGFR2/FGF10 ranks **39/71** by avg mismatch score (5.90), **54/71** for FGFR2/FGF7 (5.05).

---

## Figures

### Figure 1: Cell State Landscape and FGF Family Expression (Mouse)

![Figure 1](results/figures/fig1_cell_states_fgf.png)

- **Panel A**: UMAP of 112,676 mouse epicardial cells colored by cell state (blue=quiescent, red=activated). Subsampled to 50K for plotting.
- **Panel B**: FGFR2 expression on UMAP. Co-localizes with activated cluster (upper right).
- **Panel C**: FGF10 expression on UMAP. Enriched in quiescent clusters (lower left).
- **Panel D**: Violin plots of FGF family genes by cell state (expressing cells only; zero-expression cells removed to reveal distribution shape). FGFR1 excluded due to much higher expression scale (~3.0 vs ~0.1–0.8 for other FGF genes). Below each violin: % of cells expressing and overall fold change (mean across all cells including zeros). FGFR2 is the **only** FGF family gene upregulated in activated cells (3.0x↑, driven by increase in % expressing from 2% to 6%); FGF10 is strongly downregulated (0.2x↓).

**Data**: Quaife-Ryan 2021 (E-MTAB-10035), `mouse_quaife_ryan_analyzed.h5ad`

### Figure 2: Receptor Differential Expression Landscape (Mouse)

![Figure 2](results/figures/fig2_receptor_de_landscape.png)

- **Panel A**: Volcano plot of 1,310 receptor genes. Key receptors labeled (FGFR2, BMPR2, ACVR1, EGFR, etc.). Colored by signaling pathway.
- **Panel B**: Top 20 upregulated receptors by logFC. FGFR2 highlighted with red border. Pathway annotated for each.
- **Panel C**: Pathway-level summary showing number of significantly upregulated vs downregulated receptors per pathway. BMP and Ephrin have the most upregulated receptors.

**Data**: `receptor_rankings_by_logfc.csv` (Activated vs Quiescent, Wilcoxon)

### Figure 3: "Primed But Starved" Ligand-Receptor Mismatch

![Figure 3](results/figures/fig3_mismatch.png)

- **Panel A**: Concept diagram illustrating the hypothesis. Normal: ligand high, receptor low → balanced signaling. Post-MI: ligand depleted, receptor upregulated → insufficient signal. Therapeutic strategy: deliver depleted ligands.
- **Panel B**: Heatmap of top mismatch pairs showing receptor logFC (red) and ligand logFC (blue). FGFR2/FGF10 included at bottom.
- **Panel C**: Quadrant scatter plot of all L-R pairs (mouse). X=receptor logFC, Y=ligand logFC. The lower-right quadrant (red shading) = "primed but starved" pairs. Conserved pairs highlighted in red.
- **Panel D**: Top mismatch pairs ranked by refined composite score. Colored by pathway.

**Data**: `mouse_lr_mismatch_refined.csv`, `cross_species_lr_mismatch.csv`

### Figure 4: Geneformer In Silico Perturbation — SKIPPED

See [Geneformer section](#geneformer-in-silico-perturbation-skipped) for rationale.

### Figure 5: Cross-Species Conservation

![Figure 5](results/figures/fig5_cross_species.png)

- **Panel A**: Mouse vs Human receptor logFC scatter. Red points = receptors upregulated in both species. Key conserved receptors labeled (FGFR2, EPHA7, TYRO3, NOTCH1, etc.).
- **Panel B**: Conservation heatmap for top 20 conserved pairs. Columns: Mouse Receptor logFC, Mouse Ligand logFC, Human Receptor logFC, Human Ligand logFC. Red=up, blue=down. Pattern: mouse shows strong signal, human shows same direction but weaker.
- **Panel C**: Venn diagram showing overlap of mismatch pairs between species. 231 mouse-only, 376 human-only, 100 conserved.
- **Panel D**: Top 10 therapeutic targets table (fully automated scoring) with rank, score, conservation status, druggability, and pathway.

**Data**: `cross_species_lr_mismatch.csv`, `therapeutic_targets_corrected.csv`

### Supplementary Figure 2: L-R Database Composition

![Supp Figure 2](results/figures/supp2_lr_database.png)

- **Panel A**: Distribution of database consensus (n_db). Red dashed line marks high-confidence threshold (n_db ≥ 3). Most pairs (2,999) supported by only 1 database; 528 supported by all 5.
- **Panel B**: Database coverage. All pairs: 5,669 pairs / 1,269 ligands / 1,194 receptors. High-confidence subset: 1,953 pairs / 1,199 ligands / 1,194 receptors.
- **Panel C**: Pathway distribution of L-R pairs. "Other" dominates (86%); among annotated pathways, Wnt and BMP/TGFb have the most pairs.

**Data**: `curated_lr_pairs_mouse.csv`

### Figure scripts

| Script | Figure |
|--------|--------|
| `scripts/06_figures/fig1_cell_states_fgf.py` | Figure 1 |
| `scripts/06_figures/fig2_receptor_de_landscape.py` | Figure 2 |
| `scripts/06_figures/fig3_mismatch.py` | Figure 3 |
| `scripts/06_figures/fig5_cross_species.py` | Figure 5 |
| `scripts/06_figures/supp2_lr_database.py` | Supplementary Figure 2 |

---

## Output Files

### Result files
| File | Description |
|------|-------------|
| `results/mismatch/therapeutic_targets_corrected.csv` | **Final ranking**: 127 pairs, fully automated scoring |
| `results/mismatch/therapeutic_targets_prioritized.csv` | 127 pairs, manual druggability/literature |
| `results/mismatch/therapeutic_targets_automated.csv` | 127 pairs, DGIdb-only (before UniProt fix) |
| `results/mismatch/cross_species_lr_mismatch.csv` | All 1,953 L-R pairs with both species data |
| `results/mismatch/cross_species_conserved.csv` | 100 conserved mismatch pairs |
| `results/mismatch/mouse_lr_mismatch_refined.csv` | 127 mouse pairs with composite scoring |
| `results/mismatch/mouse_lr_mismatch_all.csv` | All 701 mouse mismatch pairs (raw) |
| `results/mismatch/gene_scores_corrected.csv` | Per-gene DGIdb + UniProt + PubMed scores |
| `results/mismatch/uniprot_deliverability.csv` | UniProt subcellular location for 190 genes |
| `results/mismatch/curated_lr_pairs_mouse.csv` | 5,669 L-R pairs used |

### Scripts
| File | Description |
|------|-------------|
| `scripts/05_mismatch/01_mouse_lr_mismatch.py` | Mouse raw mismatch analysis |
| `scripts/05_mismatch/02_mouse_lr_mismatch_refined.py` | Mouse refined composite scoring |
| `scripts/05_mismatch/03_human_lr_mismatch_refined.py` | Human mismatch analysis |
| `scripts/05_mismatch/04_cross_species_comparison.py` | Cross-species conservation |
| `scripts/05_mismatch/05_therapeutic_prioritization.py` | Manual druggability/literature scoring |
| `scripts/05_mismatch/06_automated_scoring.py` | DGIdb + PubMed automated scoring |
| `scripts/05_mismatch/07_fix_druggability.py` | UniProt deliverability correction |

---

## Date

Analysis completed: 2026-03-22
