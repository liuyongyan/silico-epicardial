# Ligand-Receptor Mismatch Analysis Results

## Objective

Systematically identify "Primed But Starved" ligand-receptor pairs in mouse epicardial cells post-MI: **receptor↑ in activated cells + cognate ligand↓**.

Positive control: FGF10/FGFR2 (wet lab validated).

---

## Data Source

- **Dataset**: Quaife-Ryan et al. 2021 (eLife), E-MTAB-10035
- **DEG comparison**: Activated vs Quiescent epicardial cells (Wilcoxon rank-sum)
- **Receptor list**: 623 significantly upregulated receptors (padj < 0.05, logFC > 0) out of 1,310 total
- **Full DEG table**: 23,376 genes with valid logFC

---

## L-R Database

**Source**: OmniPath (human, converted to mouse gene symbols via case conversion)

**Databases merged** (matching instruction 3 spec):
- CellPhoneDB
- CellTalkDB
- Ramilowski2015
- Fantom5_LRdb
- connectomeDB2020

**After complex expansion**: 5,669 unique L-R pairs (1,269 ligands × 1,194 receptors)

**Database consensus** (`n_db`): number of the 5 databases supporting each pair (range 1–5)

---

## Analysis Evolution

We tested 4 scoring approaches, progressively refining to address biases:

### Round 1: Raw Mismatch Score

`mismatch_score = receptor_logFC - ligand_logFC`

| Pair | Rank (/701) | Score |
|------|:-----------:|------:|
| Itgb1/Spon2 | 1 | 543.01 |
| Bmpr2/Bmp4 | 100 | 88.66 |
| Fgfr2/Fgf10 | **456** | 9.91 |

**Problem**: Dominated by genes with extreme logFC artifacts (Itgb1 logFC=334, biologically impossible as a true log2FC — artifact of near-zero expression in one group in scanpy's Wilcoxon test).

### Round 2: Capped logFC (±10)

Clamp logFC to biologically reasonable range.

| Pair | Rank (/701) | Score |
|------|:-----------:|------:|
| Bmpr2/Bmp4 | 82 | 20.0 |
| Fgfr2/Fgf10 | **457** | 9.91 |

**Problem**: Many pairs tied at max score (20.0). Still doesn't distinguish biologically meaningful pairs.

### Round 3: Starvation Ratio Weighting

`adjusted = raw_mismatch × starvation_ratio`

Where `starvation_ratio` = fraction of cognate ligands that are downregulated.

Key insight — the "false positives" aren't actually starved:

| Receptor | Starvation Ratio | Ligands Down/Total | Actually Starved? |
|----------|:----------------:|:------------------:|:-----------------:|
| Itgb1 | 0.394 | 39/99 | **No** — 60% of ligands upregulated, mean logFC=+40 |
| Cd44 | 0.450 | 9/20 | **No** — half of ligands upregulated |
| Sdc1 | 0.273 | 3/11 | **No** — 73% of ligands upregulated |
| **Fgfr2** | **0.652** | **15/23** | **Yes** — 65% of ligands downregulated |
| **Bmpr2** | **0.650** | **13/20** | **Yes** |
| **Acvr1** | **0.733** | **11/15** | **Yes** |

| Pair | Rank (/701) | Adjusted |
|------|:-----------:|--------:|
| Acvr1/Bmp6 | 22 | 14.67 |
| Bmpr2/Bmp4 | 34 | 13.00 |
| Fgfr2/Fgf10 | **352** | 6.46 |

**Problem**: Starvation ratio computed over ALL cognate ligands, including low-confidence pairs (n_db=1) that are often non-canonical interactions (e.g., Timp1, CCN2 for FGFR2), diluting the ratio.

### Round 4: Multi-Dimensional Scoring (Final)

```
composite = (receptor_logFC_norm) × (|ligand_logFC_norm|) × starvation_ratio × db_weight
```

**Key refinements**:
1. **logFC capped to [0, 10]** and normalized to [0, 1]
2. **Ligand must be significantly downregulated** (padj < 0.05)
3. **Starvation ratio computed using only high-confidence ligands** (n_db ≥ 3)
4. **Database consensus weight**: `db_weight = n_db / 5`

**Result**: 127 pairs (filtered from 701)

---

## Final Rankings (Round 4)

### Top 20 Pairs

| Rank | Receptor | Ligand | Composite | R logFC | L logFC | Starvation | n_db | Pathway |
|------|----------|--------|----------:|--------:|--------:|-----------:|-----:|---------|
| 1 | Unc5b | Ntn1 | 1.0000 | 59.50 | -35.44 | 1.000 | 5 | Netrin |
| 2 | Neo1 | Rgmb | 0.8000 | 49.09 | -66.80 | 1.000 | 4 | Netrin |
| 3 | Neo1 | Rgma | 0.8000 | 49.09 | -33.87 | 1.000 | 4 | Netrin |
| 4 | Neo1 | Ntn1 | 0.8000 | 49.09 | -35.44 | 1.000 | 4 | Netrin |
| 5 | Ogfr | Penk | 0.8000 | 16.17 | -129.50 | 1.000 | 4 | Opioid |
| 6 | Nectin2 | Afdn | 0.5080 | 11.86 | -6.35 | 1.000 | 4 | Adhesion |
| 7 | Itgb5 | Ccn1 | 0.4800 | 150.65 | -156.37 | 0.600 | 4 | Integrin |
| 8 | Itgb5 | Ltbp1 | 0.4800 | 150.65 | -25.43 | 0.600 | 4 | Integrin |
| 9 | Adora2b | Ntn1 | 0.4531 | 5.66 | -35.44 | 1.000 | 4 | Adenosine |
| 10 | Tyro3 | Gas6 | 0.4157 | 8.31 | -14.00 | 0.500 | 5 | TAM |
| 11 | Itgb8 | Vtn | 0.4000 | 12.00 | -120.41 | 0.500 | 4 | Integrin |
| 12 | Acvr1 | Inhbb | 0.3556 | 21.96 | -15.34 | 0.444 | 4 | BMP/Activin |
| 13 | Acvr1 | Bmp6 | 0.3556 | 21.96 | -20.73 | 0.444 | 4 | BMP/Activin |
| 14 | Plxnb2 | Sema4c | 0.3333 | 18.90 | -15.07 | 0.333 | 5 | Semaphorin |
| 15 | Bmpr2 | Bmp2 | 0.3200 | 15.12 | -12.80 | 0.400 | 4 | BMP |
| 16 | Bmpr2 | Bmp6 | 0.3200 | 15.12 | -20.73 | 0.400 | 4 | BMP |
| 17 | Bmpr2 | Bmp4 | 0.3200 | 15.12 | -73.53 | 0.400 | 4 | BMP |
| 18 | Fgfrl1 | Fgf2 | 0.3100 | 43.27 | -3.87 | 1.000 | 4 | FGF |
| 19 | Bmpr1a | Bmp4 | 0.3060 | 9.56 | -73.53 | 0.400 | 4 | BMP |
| 20 | Bmpr1a | Bmp6 | 0.3060 | 9.56 | -20.73 | 0.400 | 4 | BMP |

### Instruction 3 Key Pairs

| Pair | Inst3 Rank | **Actual Rank (/127)** | Composite | Starvation | n_db | logFC Match? |
|------|:----------:|:----------------------:|----------:|-----------:|-----:|:----:|
| Acvr1/Bmp6 | 3 | **13** | 0.3556 | 0.444 | 4 | ✓ |
| Bmpr2/Bmp4 | 1 | **17** | 0.3200 | 0.400 | 4 | ✓ |
| Bmpr1a/Bmp4 | 2 | **19** | 0.3060 | 0.400 | 4 | ✓ |
| **Fgfr2/Fgf1** | 8 | **58** | 0.1171 | 0.333 | 5 | ✓ |
| **Fgfr2/Fgf10** | 6 | **77** | 0.0817 | 0.333 | 5 | ✓ |
| **Fgfr2/Fgf7** | 7 | **78** | 0.0781 | 0.333 | 5 | ✓ |
| Fzd2/Wnt9a | 4 | filtered | — | — | <3 | — |
| Lrp6/Wnt5b | 5 | filtered | — | — | <3 | — |

**Note**: Wnt9a/Fzd2 and Wnt5b/Lrp6 were filtered out because these specific pairs have < 3 database support in our merged L-R resource.

### Best Pair per Signaling Pathway

| Pathway | Rank | Receptor | Ligand | Composite | Starvation | n_db |
|---------|:----:|----------|--------|----------:|-----------:|-----:|
| BMP/Activin | 12 | Acvr1 | Inhbb | 0.3556 | 0.444 | 4 |
| BMP | 15 | Bmpr2 | Bmp2 | 0.3200 | 0.400 | 4 |
| FGF | 18 | Fgfrl1 | Fgf2 | 0.3100 | 1.000 | 4 |
| Semaphorin | 22 | Plxna4 | Sema6a | 0.2965 | 0.500 | 4 |
| Wnt | 34 | Lrp6 | Rspo1 | 0.2000 | 0.250 | 4 |
| Sema/VEGF | 36 | Nrp2 | Pgf | 0.1925 | 0.200 | 5 |
| EGF | 51 | Egfr | Efemp1 | 0.1536 | 0.211 | 4 |
| TGF-β | 72 | Tgfbr1 | Tgfb1 | 0.0897 | 0.200 | 4 |
| Notch | 86 | Notch3 | Psen1 | 0.0547 | 0.125 | 4 |
| IGF | 113 | Igf2r | Igf2 | 0.0137 | 0.667 | 5 |

---

## Why FGFR2/FGF10 Doesn't Rank Higher

FGFR2/FGF10 ranks 77/127 — top 61% but not top 10. Three factors:

1. **Moderate receptor logFC**: FGFR2 logFC = 4.81 (rank 319/1310 among receptors). Much lower than Neo1 (49), Acvr1 (22), Bmpr2 (15), or Fzd2 (61).

2. **Moderate starvation ratio**: 5/15 high-confidence cognate ligands significantly downregulated (0.333). FGFR2 has many FGF family ligands in the database, but not all are significantly down. Compare to Neo1 (4/4 = 1.000) or Acvr1 (4/9 = 0.444).

3. **Moderate ligand logFC**: FGF10 logFC = -5.09. Compare to Bmp4 (-73.53) or Ntn1 (-35.44).

**However**, FGFR2/FGF10 has the **highest database consensus** (n_db=5, confirmed by all 5 databases), suggesting it is among the most well-established L-R pairs.

### Dimensions Not Yet Included

The following dimensions from instruction 3 Phase 6 could further improve FGFR2's ranking:
- **Geneformer perturbation effect** (Round 2 — not yet completed, requires GPU)
- **Cross-species conservation** (mouse + human)
- **Druggability** (FGF10 available as recombinant protein)
- **Literature support** (extensive evidence for FGF10/FGFR2 in cardiac repair)

---

## Comparison with Instruction 3

Instruction 3 presented FGFR2/FGF10 as "Rank 6" in its mismatch table. This was a **manually curated selection** from a few canonical signaling pathways (BMP, Wnt, FGF), not a genome-wide unbiased ranking.

Our unbiased analysis confirms:
- All logFC values in instruction 3 are **correct** (verified against our DEG CSVs)
- BMP pathway pairs (Bmpr2/Bmp4, Acvr1/Bmp6) consistently rank highest
- FGF pathway pairs rank in the middle tier by mismatch score alone
- **Simple mismatch_score alone is insufficient** to prioritize FGF10/FGFR2 — multi-dimensional scoring (Phase 6) is necessary

---

## Output Files

| File | Description |
|------|-------------|
| `results/mismatch/mouse_lr_mismatch_all.csv` | All 701 mismatch pairs (raw score) |
| `results/mismatch/mouse_lr_mismatch_top50.csv` | Top 50 by raw score |
| `results/mismatch/mouse_lr_mismatch_refined.csv` | 127 pairs with composite scoring |
| `results/mismatch/curated_lr_pairs_mouse.csv` | 5,669 L-R pairs used (mouse) |
| `scripts/05_mismatch/01_mouse_lr_mismatch.py` | Raw mismatch analysis |
| `scripts/05_mismatch/02_mouse_lr_mismatch_refined.py` | Refined multi-dimensional scoring |

---

## Date

Analysis completed: 2026-03-21
