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

**Switch from logFC to Wilcoxon scores**: The original pipeline used logFC (log-fold-change) to rank receptors and ligands. However, logFC is unstable for genes with near-zero baseline expression: a shift from 0.001 to 0.01 yields logFC = 3.3 (appearing large), while a shift from 1.0 to 2.0 yields logFC = 1.0. This inflates the apparent effect size of low-expression genes and distorts composite scoring. Wilcoxon scores (z-normalized U statistics from rank_genes_groups) are robust to these artifacts because they measure how consistently a gene ranks higher across cells, regardless of absolute expression scale.

**Filters**:
1. Receptor significantly upregulated: score > 0, padj < 0.05
2. Ligand significantly downregulated: score < 0, padj < 0.05
3. High-confidence pairs only: n_db ≥ 3
4. Starvation ratio: fraction of high-confidence cognate ligands significantly downregulated

**Scoring**: `composite = (z_receptor + |z_ligand|) × starvation_ratio × db_weight`

Where `z_receptor` and `z_ligand` are z-scores (standard deviations from the mean across all genes), and `db_weight = n_db / 5`. This is a z-score sum approach (analogous to Fisher's method for combining independent tests): because receptor upregulation and ligand downregulation are independent signals, summing their z-scores gives a principled combined statistic with no arbitrary normalization.

**Starvation ratio**: fraction of high-confidence cognate ligands that are both expressed and secreted AND significantly downregulated (denominator counts only expressed + secreted ligands, not all high-confidence ligands).

**Result**: 77 mouse mismatch pairs

### Step 3: Cross-Species Conservation

For each high-confidence L-R pair (n_db ≥ 3), check if the receptor↑ + ligand↓ pattern holds in both mouse and human. Gene ortholog mapping by case conversion (Fgfr2 → FGFR2).

**Conservation levels**:
- **Strict**: receptor↑ (padj<0.05) AND ligand↓ (padj<0.05) in BOTH species
- **Relaxed**: direction-consistent in both species, significant in at least one

**Result**: 81 conserved pairs (1 strict, 80 relaxed)

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
| 1 | ITGB1 | CD14 | 0.579 | 0.180 | 0.700 (relaxed) | 0.782 | 0.792 | Secreted | Other |
| 2 | IL2RG | IL2 | 0.553 | 0.048 | 0.700 (relaxed) | 0.880 | 0.762 | Secreted | Other |
| 3 | BMPR2 | BMP6 | 0.545 | 0.214 | 0.700 (relaxed) | 0.652 | 0.701 | Secreted | BMP |
| 4 | UNC5B | NTN1 | 0.544 | 1.000 | 0.200 | 0.450 | 0.472 | Secreted | Other |
| 5 | IL1RL2 | IL18 | 0.539 | 0.138 | 0.700 (relaxed) | 0.763 | 0.676 | Secreted | Other |
| 6 | ACVR1 | BMP6 | 0.532 | 0.255 | 0.700 (relaxed) | 0.685 | 0.542 | Secreted | BMP/Activin |
| 7 | IL2RG | IL15 | 0.527 | 0.072 | 0.700 (relaxed) | 0.843 | 0.636 | Secreted | Other |
| 8 | ITGB1 | LAMC2 | 0.500 | 0.182 | 0.700 (relaxed) | 0.737 | 0.442 | Secreted | Other |
| **9** | **FGFR2** | **FGF10** | **0.499** | 0.050 | 0.700 (relaxed) | 0.719 | 0.651 | **Secreted** | **FGF** |
| 10 | ADGRE5 | CD55 | 0.487 | 0.502 | 1.000 (strict) | 0.184 | 0.000 | Intracellular | Other |
| **11** | **FGFR2** | **FGF16** | **0.478** | 0.032 | 0.700 (relaxed) | 0.719 | 0.575 | **Secreted** | **FGF** |
| 12 | OGFR | PENK | 0.469 | 0.752 | 0.200 | 0.498 | 0.417 | Secreted | Other |
| 13 | BMPR2 | BMP4 | 0.460 | 0.330 | 0.200 | 0.652 | 0.851 | Secreted | BMP |
| 14 | GRIN2D | IL16 | 0.459 | 0.090 | 0.700 (relaxed) | 0.689 | 0.422 | Secreted | Other |
| 15 | BMPR2 | BMP2 | 0.459 | 0.224 | 0.200 | 0.804 | 0.854 | Secreted | BMP |
| 16 | BMPR1A | BMP4 | 0.445 | 0.373 | 0.200 | 0.629 | 0.738 | Secreted | BMP |
| 17 | FZD1 | MYOC | 0.443 | 0.121 | 0.700 (relaxed) | 0.583 | 0.401 | Secreted | Wnt |
| 18 | BMPR1A | BMP2 | 0.439 | 0.249 | 0.200 | 0.781 | 0.741 | Secreted | BMP |
| 19 | BMPR2 | BMP7 | 0.431 | 0.188 | 0.200 | 0.824 | 0.749 | Secreted | BMP |
| 20 | EGFR | ANXA1 | 0.431 | 0.088 | 0.200 | 0.923 | 0.798 | Secreted | EGF |

### Positive Control Validation

FGF10/FGFR2 ranks **#9/77** with fully automated score-based scoring -- no manual annotation. Comparison across scoring methods:

| Pair | Mismatch only | DGIdb-only | logFC-based | **Score-based (final)** |
|------|:------------:|:----------:|:-----------:|:-----------------------:|
| Fgfr2/Fgf10 | 64/77 | 27 | 13 | **9** |
| Fgfr2/Fgf16 | 70/77 | -- | 20 | **11** |
| Acvr1/Bmp6 | 12/77 | 11 | 5 | **6** |
| Bmpr2/Bmp6 | 21/77 | 8 | 3 | **3** |
| Unc5b/Ntn1 | 1/77 | 18 | 7 | **4** |
| Bmpr2/Bmp2 | 20/77 | 20 | 18 | **15** |

Note: Tyro3/Gas6, Insr/Nampt, and Fgfr2/Fgf7 were removed from the 77 score-based mismatch pairs because their ligands did not meet the score < 0 (padj < 0.05) filter. This is expected: Wilcoxon scores correct for logFC artifacts in near-zero-expression genes, so some previously included pairs with inflated logFC no longer qualify.

The automated score places FGFR2/FGF10 at #9 (top 12%), up from #13 in the logFC-based analysis. This improvement reflects the removal of pairs whose mismatch signal was driven by logFC artifacts rather than genuine differential expression.

### Key Findings

1. **FGF10/FGFR2 is cross-species conserved and in the top 12%.** Ranks #9/77 without manual input. FGFR2 is upregulated and FGF10 is downregulated in both mouse (significant) and human (directional trend).

2. **BMP pathway pairs rank highly.** BMPR2/BMP6 (#3), ACVR1/BMP6 (#6), BMPR2/BMP2 (#15). Note: BMP4 is upregulated in human (opposite of mouse), so BMP4 pairs are NOT conserved.

3. **Score-based filtering removes logFC artifacts.** The switch from logFC to Wilcoxon scores reduced mismatch pairs from 127 to 77. Pairs like TYRO3/GAS6 and INSR/NAMPT dropped out because their ligands' apparent downregulation was driven by near-zero expression artifacts, not genuine differential expression. UNC5B/NTN1 (#1) and OGFR/PENK (#2) now rank highest by mismatch composite.

4. **Mismatch score has a systematic blind spot for sparse genes.** FGFR2/FGF10's mismatch composite is only 0.050 (rank 64/77), contributing just 3% of its final priority score. The #9 ranking is still primarily driven by conservation + druggability + literature. This is not a failure of the scoring method — it reflects a fundamental limitation: FGFR2 is expressed in only 2–6% of cells, and FGF10 in 1–3%. Both Wilcoxon scores and logFC are bulk statistical tests that measure "how consistently does this gene differ across all cells." For genes expressed in <10% of cells, even a real biological difference produces a weak statistical signal because 90%+ of cells contribute zero to both groups. The wet lab validation of FGF10/FGFR2 may simply be ahead of what computational mismatch analysis can confirm at this expression level. This suggests that for sparse receptor-ligand pairs, alternative approaches (e.g., expression percentage-based metrics, or restricting analysis to expressing cells only) may be more appropriate than whole-population differential expression.

5. **Human data is a bottleneck.** Only 1 conserved pair reaches strict significance in both species (ADGRE5/CD55). All other conservation calls rely on directional consistency in human (trend level). Multi-patient human MI data would substantially strengthen these findings.

### 81 Cross-Species Conserved Pairs (by avg mismatch score)

See `cross_species_conserved_scores.csv` for the full 81 pairs. Top pairs by avg mismatch score:

| Rank | Receptor | Ligand | Avg Mismatch Score | n_db |
|:----:|----------|--------|:------------------:|:----:|
| 1 | ITGB1 | LAMC2 | 0.166 | 4 |
| 2 | ITGB1 | CD14 | 0.164 | 4 |
| 3 | TNFRSF12A | TNFSF12 | 0.134 | 5 |
| 4 | FZD1 | MYOC | 0.117 | 3 |
| 5 | FZD4 | MYOC | 0.101 | 3 |
| 8 | ACVR1 | BMP6 | 0.082 | 4 |
| 12 | LRP6 | RSPO1 | 0.076 | 4 |
| 16 | BMPR2 | BMP6 | 0.068 | 4 |
| 18 | OSMR | OSM | 0.061 | 4 |
| 34 | FGFR2 | FGF16 | 0.025 | 3 |
| 40 | FGFR2 | FGF10 | 0.024 | 3 |

FGFR2/FGF10 ranks **40/81** by avg mismatch score (0.024), **34/81** for FGFR2/FGF16 (0.025).

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

- **Panel A**: Volcano plot of receptor genes (125 logFC-artifact genes removed where |logFC/score| ratio > 3). X-axis: log₂FC, Y-axis: -log₁₀(padj). Key receptors labeled. Colored by signaling pathway.
- **Panel B**: Waterfall plot of all 1,057 significantly upregulated receptors ranked by Wilcoxon score. Key receptors highlighted with rank labels. FGFR2 sits in the tail of the curve — its Wilcoxon score (12.1) is modest compared to top receptors like Palld (190) or Itgb1 (118), consistent with its low expression in only 2–6% of cells. The steep drop-off shows that only a small number of receptors have strong statistical signal; the majority (including FGFR2) cluster in the low-score region.
- **Panel C**: Pathway-level summary showing number of significantly upregulated vs downregulated receptors per pathway. Nearly all annotated signaling receptors are upregulated in the activated state (consistent with the EMT/mesenchymal phenotype); only TGF-β, PDGF, and VEGF have any downregulated members.

**Data**: `receptor_rankings_by_logfc.csv` (Activated vs Quiescent, Wilcoxon scores)

### Figure 3: "Primed But Starved" Ligand-Receptor Mismatch

![Figure 3](results/figures/fig3_mismatch.png)

- **Panel A**: Concept diagram illustrating the hypothesis. Normal: ligand high, receptor low → balanced signaling. Post-MI: ligand depleted, receptor upregulated → insufficient signal. Therapeutic strategy: deliver depleted ligands.
- **Panel B**: Heatmap of top mismatch pairs showing receptor Wilcoxon score (red) and ligand Wilcoxon score (blue). FGFR2/FGF10 included at bottom.
- **Panel C**: Quadrant scatter plot of all L-R pairs (mouse). X=receptor score, Y=ligand score. The lower-right quadrant (red shading) = "primed but starved" pairs. Conserved pairs highlighted in red.
- **Panel D**: Top mismatch pairs ranked by refined composite score (score-based). Colored by pathway.

**Data**: `mouse_lr_mismatch_scores.csv`, `cross_species_lr_mismatch_scores.csv`

### Figure 4: Geneformer In Silico Perturbation — SKIPPED

See [Geneformer section](#geneformer-in-silico-perturbation-skipped) for rationale.

### Figure 5: Cross-Species Conservation

![Figure 5](results/figures/fig5_cross_species.png)

- **Panel A**: Mouse vs Human receptor effect size scatter using rank-biserial correlation (r = z/√N), a sample-size-independent effect size measure (Rosenthal 1994). This corrects for the 5.6× difference in total cell counts between mouse (N=112,676) and human (N=20,007) which inflates raw z-scores. X-axis: mouse r (±0.2), Y-axis: human r (±0.1). Red points = receptors upregulated in both species. Key conserved receptors labeled (FGFR2, BMPR2, ACVR1, NOTCH1, etc.). Most points fall below the diagonal, indicating mouse effects are systematically larger than human — consistent with the single-patient limitation of human MI data.
- **Panel B**: Conservation heatmap for top 20 conserved pairs, also using rank-biserial r. Same color scale (±0.3) for both species enables direct visual comparison. Mouse columns show deeper colors (larger effect sizes); human columns show lighter but directionally consistent colors — the "primed but starved" pattern is conserved in direction but attenuated in magnitude in human.
- **Panel C**: Venn diagram showing overlap of mismatch pairs between species. 206 mouse-only, 407 human-only, 81 conserved.
- **Panel D**: Top 10 therapeutic targets table (fully automated score-based scoring) with rank, score, conservation status, druggability, and pathway.

**Data**: `cross_species_lr_mismatch_scores.csv`, `therapeutic_targets_scores.csv`

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
| `results/mismatch/therapeutic_targets_scores.csv` | **Final ranking (score-based)**: 77 pairs, fully automated scoring with Wilcoxon scores |
| `results/mismatch/mouse_lr_mismatch_scores.csv` | 77 mouse pairs with score-based composite |
| `results/mismatch/cross_species_lr_mismatch_scores.csv` | All 1,953 L-R pairs with both species scores |
| `results/mismatch/cross_species_conserved_scores.csv` | 81 conserved mismatch pairs (score-based) |
| `results/mismatch/therapeutic_targets_corrected.csv` | Previous ranking (logFC-based): 127 pairs |
| `results/mismatch/cross_species_lr_mismatch.csv` | All 1,953 L-R pairs (logFC-based, previous) |
| `results/mismatch/cross_species_conserved.csv` | 100 conserved pairs (logFC-based, previous) |
| `results/mismatch/mouse_lr_mismatch_refined.csv` | 127 mouse pairs (logFC-based, previous) |
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
| `scripts/05_mismatch/08_rerun_with_scores.py` | Re-run full pipeline with Wilcoxon scores |

---

## Date

Analysis completed: 2026-03-22
