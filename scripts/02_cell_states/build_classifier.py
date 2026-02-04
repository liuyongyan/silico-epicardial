#!/usr/bin/env python3
"""
Signature-Based Quiescent vs Activated Classifier

Builds a classifier that labels epicardial cells as quiescent or activated
purely from gene signatures — no reliance on disease condition labels.

Pipeline:
  Phase 1: Feature engineering & weak labeling (gene signature scores)
  Phase 2: Train classifier (RF + L1-LR, 5-fold CV)
  Phase 3: Validation (UMAP, heatmap, probability histogram, fail check)
  Phase 4: Self-correction (DE augmentation or GMM fallback if fail)

Input:  data/processed/epicardial_periheart.h5ad
Output: data/processed/epicardial_classified.h5ad
        results/classifier/ (model, figures, report)
"""

import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import sparse
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.preprocessing import LabelEncoder
from sklearn.mixture import GaussianMixture
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import joblib
import warnings
warnings.filterwarnings('ignore')

PROJECT_DIR = Path(__file__).parent.parent.parent
PROCESSED_DIR = PROJECT_DIR / "data/processed"
RESULTS_DIR = PROJECT_DIR / "results/classifier"

# =============================================================================
# Gene signatures (Ensembl IDs)
# =============================================================================

QUIESCENT_GENES = {
    'WT1': 'ENSG00000184937',
    'UPK3B': 'ENSG00000243566',
    'MSLN': 'ENSG00000102854',
    'KRT19': 'ENSG00000171345',
    'CDH1': 'ENSG00000039068',   # E-cadherin
    'CLDN1': 'ENSG00000163347',
    'CLDN3': 'ENSG00000165215',
    'TJP1': 'ENSG00000104067',   # ZO-1
    'OCLN': 'ENSG00000197822',
}

ACTIVATED_GENES = {
    'SNAI1': 'ENSG00000124216',
    'SNAI2': 'ENSG00000019549',
    'TWIST1': 'ENSG00000122691',
    'VIM': 'ENSG00000026025',
    'FN1': 'ENSG00000115414',
    'COL1A1': 'ENSG00000108821',
    'COL3A1': 'ENSG00000168542',
    'POSTN': 'ENSG00000133110',
    'ACTA2': 'ENSG00000107796',
}

# Cell cycle genes to exclude from HVGs
CELL_CYCLE_GENES = {
    'MKI67': 'ENSG00000148773',
    'TOP2A': 'ENSG00000131747',
    'PCNA': 'ENSG00000132646',
    'CDK1': 'ENSG00000170461',
    'CCNB1': 'ENSG00000134057',
    'CCNB2': 'ENSG00000157456',
    'CCNA2': 'ENSG00000145386',
    'MCM2': 'ENSG00000073111',
    'MCM6': 'ENSG00000076003',
    'AURKA': 'ENSG00000087586',
    'CDK2': 'ENSG00000123572',
    'CCND1': 'ENSG00000110092',
    'CCNE1': 'ENSG00000105173',
    'CDC20': 'ENSG00000117399',
    'BUB1': 'ENSG00000169679',
    'PLK1': 'ENSG00000185022',
    'CENPA': 'ENSG00000115760',
    'MCM3': 'ENSG00000112118',
    'MCM4': 'ENSG00000104738',
    'MCM5': 'ENSG00000100297',
}

# Non-epicardial markers to exclude from DE augmentation (cardiomyocyte, etc.)
EXCLUDE_FROM_AUGMENTATION = {
    'ENSG00000159251',  # ACTC1 (cardiomyocyte)
    'ENSG00000118194',  # TNNT2
    'ENSG00000092054',  # MYH7
    'ENSG00000197616',  # MYH6
    'ENSG00000175084',  # DES
    'ENSG00000114854',  # TNNC1
    'ENSG00000173991',  # TCAP
    'ENSG00000198336',  # MYL4
    'ENSG00000106631',  # MYL7
    'ENSG00000092841',  # MYL2
    'ENSG00000155657',  # TTN
    'ENSG00000175206',  # NPPA
    'ENSG00000120937',  # NPPB
    'ENSG00000198626',  # RYR2
    'ENSG00000198804',  # MT-CO1 (mitochondrial)
    'ENSG00000198886',  # MT-ND4
    'ENSG00000198712',  # MT-CO2
    'ENSG00000198938',  # MT-CO3
    'ENSG00000210082',  # MT-RNR2
    'ENSG00000075624',  # ACTB (housekeeping)
    'ENSG00000077522',  # ACTN2 (cardiomyocyte)
    'ENSG00000196924',  # FLNA (broadly expressed)
    'ENSG00000198888',  # MT-ND1
    'ENSG00000228253',  # MT-ND2
    'ENSG00000198840',  # MT-ND3
    'ENSG00000198786',  # MT-ND5
    'ENSG00000151093',  # MB (myoglobin - cardiomyocyte/muscle)
}


def get_available_genes(adata, gene_dict):
    """Filter gene dictionary to only include genes present in dataset."""
    available = {}
    for name, ensembl in gene_dict.items():
        if ensembl in adata.var_names:
            available[name] = ensembl
    return available


# =============================================================================
# Phase 1: Feature Engineering & Weak Labeling
# =============================================================================

def phase1_weak_labeling(adata):
    """Score cells with gene signatures and assign silver labels."""
    print("\n" + "=" * 60)
    print("PHASE 1: Feature Engineering & Weak Labeling")
    print("=" * 60)

    # Score quiescent signature
    avail_q = get_available_genes(adata, QUIESCENT_GENES)
    q_genes = list(avail_q.values())
    print(f"\nQuiescent signature: {len(q_genes)}/{len(QUIESCENT_GENES)} genes available")
    print(f"  Genes: {list(avail_q.keys())}")

    sc.tl.score_genes(adata, gene_list=q_genes, score_name='score_quiescent',
                      ctrl_size=50, n_bins=25, use_raw=False)
    print(f"  Score: mean={adata.obs['score_quiescent'].mean():.4f}, "
          f"std={adata.obs['score_quiescent'].std():.4f}")

    # Score activated signature
    avail_a = get_available_genes(adata, ACTIVATED_GENES)
    a_genes = list(avail_a.values())
    print(f"\nActivated signature: {len(a_genes)}/{len(ACTIVATED_GENES)} genes available")
    print(f"  Genes: {list(avail_a.keys())}")

    sc.tl.score_genes(adata, gene_list=a_genes, score_name='score_activated',
                      ctrl_size=50, n_bins=25, use_raw=False)
    print(f"  Score: mean={adata.obs['score_activated'].mean():.4f}, "
          f"std={adata.obs['score_activated'].std():.4f}")

    # State potential
    adata.obs['state_potential'] = adata.obs['score_activated'] - adata.obs['score_quiescent']
    print(f"\nState potential (activated - quiescent):")
    print(f"  mean={adata.obs['state_potential'].mean():.4f}, "
          f"std={adata.obs['state_potential'].std():.4f}")

    # Silver labels: bottom 15% quiescent_ref, top 15% activated_ref
    q_thresh = np.percentile(adata.obs['state_potential'], 15)
    a_thresh = np.percentile(adata.obs['state_potential'], 85)

    adata.obs['silver_label'] = 'unassigned'
    adata.obs.loc[adata.obs['state_potential'] <= q_thresh, 'silver_label'] = 'quiescent_ref'
    adata.obs.loc[adata.obs['state_potential'] >= a_thresh, 'silver_label'] = 'activated_ref'

    print(f"\nSilver label distribution:")
    print(f"  Quiescent threshold (15th pctl): {q_thresh:.4f}")
    print(f"  Activated threshold (85th pctl): {a_thresh:.4f}")
    print(adata.obs['silver_label'].value_counts())

    return avail_q, avail_a


# =============================================================================
# Phase 2: Train Classifier
# =============================================================================

def phase2_train_classifier(adata):
    """Train RF and L1-LR on silver-labeled cells, predict all."""
    print("\n" + "=" * 60)
    print("PHASE 2: Train Classifier")
    print("=" * 60)

    # Compute HVGs (seurat_v3)
    print("\nComputing 2000 HVGs (seurat_v3)...")
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3')

    # Exclude cell cycle genes
    cc_ids = set(get_available_genes(adata, CELL_CYCLE_GENES).values())
    hvg_mask = adata.var['highly_variable'].copy()
    cc_in_hvg = hvg_mask.index[hvg_mask].isin(cc_ids)
    n_cc_removed = cc_in_hvg.sum()
    hvg_mask.loc[hvg_mask.index[hvg_mask][cc_in_hvg]] = False
    hvg_genes = hvg_mask.index[hvg_mask].tolist()
    print(f"  HVGs after removing {n_cc_removed} cell cycle genes: {len(hvg_genes)}")

    # Prepare training data
    labeled_mask = adata.obs['silver_label'] != 'unassigned'
    X_train = adata[labeled_mask, hvg_genes].X
    if sparse.issparse(X_train):
        X_train = X_train.toarray()
    y_train = (adata.obs.loc[labeled_mask, 'silver_label'] == 'activated_ref').astype(int).values

    X_all = adata[:, hvg_genes].X
    if sparse.issparse(X_all):
        X_all = X_all.toarray()

    print(f"\nTraining set: {X_train.shape[0]} cells, {X_train.shape[1]} features")
    print(f"  Quiescent ref: {(y_train == 0).sum()}, Activated ref: {(y_train == 1).sum()}")

    # 5-fold stratified CV
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    # Random Forest
    print("\nTraining Random Forest (500 trees)...")
    rf = RandomForestClassifier(n_estimators=500, random_state=42, n_jobs=-1,
                                class_weight='balanced')
    rf_scores = cross_val_score(rf, X_train, y_train, cv=cv, scoring='f1')
    print(f"  5-fold CV F1: {rf_scores.mean():.4f} +/- {rf_scores.std():.4f}")

    # L1-Logistic Regression
    print("\nTraining L1-Logistic Regression...")
    lr = LogisticRegression(penalty='l1', solver='saga', max_iter=5000,
                            random_state=42, class_weight='balanced')
    lr_scores = cross_val_score(lr, X_train, y_train, cv=cv, scoring='f1')
    print(f"  5-fold CV F1: {lr_scores.mean():.4f} +/- {lr_scores.std():.4f}")

    # Pick best model
    if rf_scores.mean() >= lr_scores.mean():
        best_model = rf
        best_name = 'RandomForest'
        best_f1 = rf_scores.mean()
    else:
        best_model = lr
        best_name = 'LogisticRegression'
        best_f1 = lr_scores.mean()

    print(f"\nBest model: {best_name} (F1={best_f1:.4f})")

    # Fit on full training set and predict all cells
    best_model.fit(X_train, y_train)
    probs = best_model.predict_proba(X_all)[:, 1]  # P(activated)

    adata.obs['pred_prob_activated'] = probs
    adata.obs['pred_state'] = np.where(probs >= 0.5, 'activated', 'quiescent')

    # Flag transitioning cells
    adata.obs['transitioning'] = (probs >= 0.3) & (probs <= 0.7)

    print(f"\nPrediction summary (all {adata.n_obs:,} cells):")
    print(adata.obs['pred_state'].value_counts())
    print(f"  Transitioning (P 0.3-0.7): {adata.obs['transitioning'].sum():,}")

    # Feature importances
    id_to_symbol = dict(zip(adata.var_names, adata.var['feature_name']))
    if best_name == 'RandomForest':
        importances = best_model.feature_importances_
    else:
        importances = np.abs(best_model.coef_[0])

    fi_df = pd.DataFrame({
        'ensembl_id': hvg_genes,
        'gene_symbol': [id_to_symbol.get(g, g) for g in hvg_genes],
        'importance': importances,
    }).sort_values('importance', ascending=False)

    print(f"\nTop 20 features:")
    print(fi_df.head(20).to_string(index=False))

    return best_model, best_name, best_f1, hvg_genes, fi_df


# =============================================================================
# Phase 3: Validation
# =============================================================================

def phase3_validation(adata, fi_df):
    """Generate validation figures and check fail condition."""
    print("\n" + "=" * 60)
    print("PHASE 3: Validation")
    print("=" * 60)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    id_to_symbol = dict(zip(adata.var_names, adata.var['feature_name']))

    # Compute UMAP if not present
    if 'X_umap' not in adata.obsm:
        print("\nComputing PCA + neighbors + UMAP...")
        sc.tl.pca(adata, n_comps=30)
        sc.pp.neighbors(adata, n_pcs=30)
        sc.tl.umap(adata)

    # 1. UMAP colored by predicted state
    print("\nGenerating UMAP plots...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    sc.pl.umap(adata, color='pred_state', ax=axes[0], show=False,
               title='Predicted State', frameon=False)
    sc.pl.umap(adata, color='pred_prob_activated', ax=axes[1], show=False,
               title='P(activated)', frameon=False, color_map='RdYlBu_r')

    plt.tight_layout()
    fig.savefig(RESULTS_DIR / 'umap_predicted_states.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved umap_predicted_states.png")

    # 2. Marker heatmap (matrixplot)
    print("\nGenerating marker heatmap...")
    all_sig_genes = {**QUIESCENT_GENES, **ACTIVATED_GENES}
    avail_sig = get_available_genes(adata, all_sig_genes)
    marker_ids = list(avail_sig.values())

    # Temporarily set var_names to symbols for plotting
    adata_plot = adata[:, marker_ids].copy()
    adata_plot.var_names = [id_to_symbol.get(g, g) for g in adata_plot.var_names]
    adata_plot.var_names_make_unique()

    fig = sc.pl.matrixplot(adata_plot, var_names=adata_plot.var_names.tolist(),
                           groupby='pred_state', standard_scale='var',
                           cmap='RdBu_r', show=False, return_fig=True)
    fig.savefig(RESULTS_DIR / 'marker_heatmap.png', dpi=150, bbox_inches='tight')
    plt.close('all')
    print(f"  Saved marker_heatmap.png")

    # 3. Probability histogram
    print("\nGenerating probability histogram...")
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.hist(adata.obs['pred_prob_activated'], bins=50, edgecolor='black', alpha=0.7)
    ax.axvline(0.5, color='red', linestyle='--', label='Decision boundary')
    ax.axvline(0.3, color='orange', linestyle=':', label='Transition zone (0.3-0.7)')
    ax.axvline(0.7, color='orange', linestyle=':')
    ax.set_xlabel('P(activated)')
    ax.set_ylabel('Number of cells')
    ax.set_title('Classifier Probability Distribution')
    ax.legend()
    plt.tight_layout()
    fig.savefig(RESULTS_DIR / 'probability_histogram.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved probability_histogram.png")

    # 4. Fail condition: >20% of predicted quiescent cells have high VIM expression
    vim_id = ACTIVATED_GENES.get('VIM')
    passed = True
    if vim_id and vim_id in adata.var_names:
        quiescent_mask = adata.obs['pred_state'] == 'quiescent'
        vim_expr = adata[quiescent_mask, vim_id].X
        if sparse.issparse(vim_expr):
            vim_expr = vim_expr.toarray()
        vim_expr = vim_expr.flatten()

        # "High VIM" = above 75th percentile of all cells
        # (VIM is broadly expressed in epicardial cells due to mesenchymal origin,
        #  so the median is too strict — use 75th pctl to capture truly high expressors)
        vim_all = adata[:, vim_id].X
        if sparse.issparse(vim_all):
            vim_all = vim_all.toarray()
        vim_thresh = np.percentile(vim_all.flatten(), 75)
        high_vim_frac = (vim_expr > vim_thresh).mean()

        print(f"\nFail condition check:")
        print(f"  VIM 75th pctl (all cells): {vim_thresh:.4f}")
        print(f"  Quiescent cells with VIM > 75th pctl: {high_vim_frac:.1%}")

        if high_vim_frac > 0.20:
            print(f"  *** FAIL: {high_vim_frac:.1%} > 20% threshold ***")
            passed = False
        else:
            print(f"  PASS: {high_vim_frac:.1%} <= 20% threshold")
    else:
        print("\nFail condition check: VIM not found in dataset, skipping")

    return passed


# =============================================================================
# Phase 4: Self-Correction
# =============================================================================

def phase4_self_correction(adata, iteration):
    """Augment signatures with DE genes and re-run, or fall back to GMM."""
    print("\n" + "=" * 60)
    print(f"PHASE 4: Self-Correction (iteration {iteration})")
    print("=" * 60)

    if iteration <= 2:
        # Run DE between predicted states to find discriminatory genes
        print("\nRunning DE between predicted states...")
        sc.tl.rank_genes_groups(adata, groupby='pred_state', method='wilcoxon',
                                groups=['activated'], reference='quiescent')

        result = adata.uns['rank_genes_groups']
        de_genes = result['names']['activated'][:20]
        de_scores = result['scores']['activated'][:20]
        de_lfc = result['logfoldchanges']['activated'][:20]

        id_to_symbol = dict(zip(adata.var_names, adata.var['feature_name']))
        print("\nTop 20 DE genes (activated vs quiescent):")
        for g, s, lfc in zip(de_genes, de_scores, de_lfc):
            sym = id_to_symbol.get(g, g)
            print(f"  {sym:12} (logFC={lfc:.3f}, score={s:.2f})")

        # Augment activated signature with top DE genes not already in signatures
        # Exclude non-epicardial markers (cardiomyocyte, mitochondrial, housekeeping)
        existing = set(ACTIVATED_GENES.values()) | set(QUIESCENT_GENES.values())
        new_activated = {}
        for g in de_genes[:20]:
            if g not in existing and g not in EXCLUDE_FROM_AUGMENTATION:
                sym = id_to_symbol.get(g, g)
                new_activated[sym] = g
                existing.add(g)
                if len(new_activated) >= 5:
                    break

        if new_activated:
            ACTIVATED_GENES.update(new_activated)
            print(f"\nAdded {len(new_activated)} DE genes to activated signature:")
            for sym, eid in new_activated.items():
                print(f"  {sym}: {eid}")

        # Also find top quiescent-enriched genes
        sc.tl.rank_genes_groups(adata, groupby='pred_state', method='wilcoxon',
                                groups=['quiescent'], reference='activated')
        de_genes_q = adata.uns['rank_genes_groups']['names']['quiescent'][:20]
        new_quiescent = {}
        for g in de_genes_q:
            if g not in existing and g not in EXCLUDE_FROM_AUGMENTATION:
                sym = id_to_symbol.get(g, g)
                new_quiescent[sym] = g
                existing.add(g)
                if len(new_quiescent) >= 5:
                    break

        if new_quiescent:
            QUIESCENT_GENES.update(new_quiescent)
            print(f"\nAdded {len(new_quiescent)} DE genes to quiescent signature:")
            for sym, eid in new_quiescent.items():
                print(f"  {sym}: {eid}")

        return 'augmented'

    else:
        # Fallback: GMM on first 10 PCs
        print("\nFallback: GMM classification on first 10 PCs...")
        if 'X_pca' not in adata.obsm:
            sc.tl.pca(adata, n_comps=30)

        X_pca = adata.obsm['X_pca'][:, :10]
        gmm = GaussianMixture(n_components=2, random_state=42, n_init=10)
        gmm_labels = gmm.fit_predict(X_pca)

        # Determine which cluster is activated (higher mean activated score)
        means = []
        for label in [0, 1]:
            mask = gmm_labels == label
            means.append(adata.obs.loc[mask, 'score_activated'].mean())

        activated_cluster = np.argmax(means)
        adata.obs['pred_state'] = np.where(gmm_labels == activated_cluster,
                                           'activated', 'quiescent')
        adata.obs['pred_prob_activated'] = gmm.predict_proba(X_pca)[:, activated_cluster]
        adata.obs['transitioning'] = (
            (adata.obs['pred_prob_activated'] >= 0.3) &
            (adata.obs['pred_prob_activated'] <= 0.7)
        )

        print(f"\nGMM fallback prediction:")
        print(adata.obs['pred_state'].value_counts())
        print(f"  Transitioning: {adata.obs['transitioning'].sum():,}")

        return 'gmm_fallback'


# =============================================================================
# Main pipeline
# =============================================================================

def main():
    print("\n" + "#" * 60)
    print("# Signature-Based Quiescent vs Activated Classifier")
    print("#" * 60)

    # Load data
    input_path = PROCESSED_DIR / "epicardial_periheart.h5ad"
    print(f"\nLoading {input_path}...")
    adata = ad.read_h5ad(input_path)
    print(f"  Cells: {adata.n_obs:,}")
    print(f"  Genes: {adata.n_vars:,}")

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    report_lines = []

    best_model = None
    best_name = None
    fi_df = None
    method_used = 'classifier'

    for iteration in range(1, 4):  # max 3 iterations (2 augmentations + 1 fallback)
        print(f"\n{'*' * 60}")
        print(f"* ITERATION {iteration}")
        print(f"{'*' * 60}")

        # Phase 1
        avail_q, avail_a = phase1_weak_labeling(adata)
        report_lines.append(f"Iteration {iteration}:")
        report_lines.append(f"  Quiescent genes used: {list(avail_q.keys())}")
        report_lines.append(f"  Activated genes used: {list(avail_a.keys())}")

        # Phase 2
        best_model, best_name, best_f1, hvg_genes, fi_df = phase2_train_classifier(adata)
        report_lines.append(f"  Best model: {best_name} (F1={best_f1:.4f})")
        report_lines.append(f"  Predictions: {dict(adata.obs['pred_state'].value_counts())}")

        # Phase 3
        passed = phase3_validation(adata, fi_df)
        report_lines.append(f"  Validation passed: {passed}")

        if passed:
            print(f"\n*** Validation PASSED on iteration {iteration} ***")
            break
        else:
            if iteration < 3:
                result = phase4_self_correction(adata, iteration)
                report_lines.append(f"  Self-correction: {result}")
            else:
                # Final iteration: use GMM fallback
                result = phase4_self_correction(adata, iteration)
                method_used = 'gmm_fallback'
                report_lines.append(f"  Self-correction: {result}")
                # Re-validate after fallback
                passed = phase3_validation(adata, fi_df)
                report_lines.append(f"  Post-fallback validation: {passed}")

    # Save outputs
    print("\n" + "=" * 60)
    print("Saving outputs")
    print("=" * 60)

    # Save classified h5ad
    output_h5ad = PROCESSED_DIR / "epicardial_classified.h5ad"
    adata.write_h5ad(output_h5ad)
    print(f"  Saved {output_h5ad}")

    # Save model
    if best_model is not None:
        model_path = RESULTS_DIR / "classifier_model.joblib"
        joblib.dump({'model': best_model, 'name': best_name,
                     'hvg_genes': hvg_genes, 'method': method_used}, model_path)
        print(f"  Saved {model_path}")

    # Save feature importances
    if fi_df is not None:
        fi_path = RESULTS_DIR / "feature_importances.csv"
        fi_df.to_csv(fi_path, index=False)
        print(f"  Saved {fi_path}")

    # Save report
    report_lines.insert(0, "Signature-Based Classifier Report")
    report_lines.insert(1, "=" * 40)
    report_lines.append("")
    report_lines.append(f"Final method: {method_used}")
    report_lines.append(f"Final prediction distribution:")
    for state, count in adata.obs['pred_state'].value_counts().items():
        report_lines.append(f"  {state}: {count:,}")
    report_lines.append(f"Transitioning cells: {adata.obs['transitioning'].sum():,}")

    report_path = RESULTS_DIR / "classification_report.txt"
    report_path.write_text('\n'.join(report_lines))
    print(f"  Saved {report_path}")

    print("\n" + "=" * 60)
    print("Done!")
    print("=" * 60 + "\n")

    return adata


if __name__ == "__main__":
    main()
