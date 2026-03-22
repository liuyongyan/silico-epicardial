#!/usr/bin/env python3
"""
Automated Druggability and Literature Scoring via DGIdb + PubMed APIs.

Druggability: DGIdb GraphQL API — number of drug interactions, approved drugs
Literature: PubMed E-utilities — publication count for gene + cardiac terms

Then re-run Phase 6 prioritization with automated scores.
"""

import pandas as pd
import numpy as np
import json
import urllib.request
import time
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent.parent.parent
MISMATCH_DIR = PROJECT_DIR / "results" / "mismatch"

# ============================================================================
# 1. Get unique genes from conserved pairs
# ============================================================================
print("=" * 80)
print("Step 1: Collecting genes to query")
print("=" * 80)

conserved = pd.read_csv(MISMATCH_DIR / "cross_species_conserved.csv")
mouse_refined = pd.read_csv(MISMATCH_DIR / "mouse_lr_mismatch_refined.csv")

# All genes from both datasets
all_genes = set()
all_genes.update(conserved['receptor'].tolist())
all_genes.update(conserved['ligand'].tolist())
all_genes.update(mouse_refined['receptor'].str.upper().tolist())
all_genes.update(mouse_refined['ligand'].str.upper().tolist())
all_genes = sorted(all_genes)
print(f"Total unique genes: {len(all_genes)}")

# ============================================================================
# 2. Query DGIdb for druggability
# ============================================================================
print("\n" + "=" * 80)
print("Step 2: Querying DGIdb for drug-gene interactions")
print("=" * 80)

def query_dgidb(gene_list, batch_size=25):
    """Query DGIdb GraphQL API in batches."""
    results = {}
    url = 'https://dgidb.org/api/graphql'

    for i in range(0, len(gene_list), batch_size):
        batch = gene_list[i:i+batch_size]
        names_str = ', '.join(f'"{g}"' for g in batch)

        query = f'''
        {{
          genes(names: [{names_str}]) {{
            nodes {{
              name
              interactions {{
                drug {{
                  name
                  approved
                }}
                interactionTypes {{
                  type
                  directionality
                }}
              }}
            }}
          }}
        }}
        '''

        try:
            data = json.dumps({'query': query}).encode()
            req = urllib.request.Request(url, data=data, headers={
                'Content-Type': 'application/json',
                'User-Agent': 'EpicardialAnalysis/1.0'
            })
            resp = urllib.request.urlopen(req, timeout=30)
            result = json.loads(resp.read())

            for gene in result['data']['genes']['nodes']:
                name = gene['name']
                interactions = gene['interactions']
                n_total = len(interactions)
                approved_drugs = [x['drug']['name'] for x in interactions
                                  if x['drug'].get('approved')]
                n_approved = len(approved_drugs)

                # Collect interaction types
                types = set()
                for inter in interactions:
                    for t in inter['interactionTypes']:
                        if t['type']:
                            types.add(t['type'])

                has_agonist = 'agonist' in types or 'activator' in types
                has_antibody = 'antibody' in types

                results[name] = {
                    'n_interactions': n_total,
                    'n_approved': n_approved,
                    'approved_drugs': '; '.join(approved_drugs[:5]),
                    'interaction_types': '; '.join(sorted(types)),
                    'has_agonist': has_agonist,
                    'has_antibody': has_antibody,
                }

            print(f"  Batch {i//batch_size+1}: queried {len(batch)} genes, "
                  f"found {len([g for g in batch if g in results])}")

        except Exception as e:
            print(f"  Batch {i//batch_size+1}: ERROR - {e}")

        time.sleep(0.5)  # rate limiting

    return results

dgidb_results = query_dgidb(all_genes)
print(f"\nDGIdb results: {len(dgidb_results)} genes with data")

# ============================================================================
# 3. Query PubMed for literature support
# ============================================================================
print("\n" + "=" * 80)
print("Step 3: Querying PubMed for cardiac literature counts")
print("=" * 80)

def query_pubmed(gene_list, context_terms=None):
    """Query PubMed for publication counts."""
    if context_terms is None:
        context_terms = '(cardiac OR epicardial OR "heart repair" OR "myocardial infarction" OR "cardiac regeneration")'

    results = {}
    base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'

    for i, gene in enumerate(gene_list):
        query = f'{gene}+AND+{context_terms}'.replace(' ', '+')
        url = f'{base_url}?db=pubmed&term={query}&retmode=json&retmax=0'

        try:
            req = urllib.request.Request(url, headers={
                'User-Agent': 'EpicardialAnalysis/1.0'
            })
            resp = urllib.request.urlopen(req, timeout=10)
            data = json.loads(resp.read())
            count = int(data['esearchresult']['count'])
            results[gene] = count
        except Exception as e:
            results[gene] = 0
            print(f"  {gene}: ERROR - {e}")

        # Rate limiting (NCBI allows 3 requests/sec without API key)
        if (i + 1) % 3 == 0:
            time.sleep(0.4)

        if (i + 1) % 25 == 0:
            print(f"  Queried {i+1}/{len(gene_list)} genes...")

    return results

pubmed_results = query_pubmed(all_genes)
print(f"\nPubMed results: {len(pubmed_results)} genes queried")

# Show top 20 by publication count
pubmed_df = pd.DataFrame([
    {'gene': k, 'pubmed_cardiac': v}
    for k, v in sorted(pubmed_results.items(), key=lambda x: -x[1])
])
print("\nTop 20 genes by cardiac literature:")
print(pubmed_df.head(20).to_string(index=False))

# ============================================================================
# 4. Build automated scores
# ============================================================================
print("\n" + "=" * 80)
print("Step 4: Computing automated druggability and literature scores")
print("=" * 80)

# Druggability score:
# Based on number of drug interactions and approved drugs
# Normalize: log-scale to handle wide range (0 to 87 interactions)
def druggability_score(gene):
    if gene not in dgidb_results:
        return 0.1  # no data
    d = dgidb_results[gene]
    # Score components:
    # - Has any interaction: +0.2
    # - log(n_interactions+1) normalized: up to +0.4
    # - Has approved drugs: +0.2
    # - Has agonist (most relevant for "deliver ligand" strategy): +0.2
    score = 0.0
    if d['n_interactions'] > 0:
        score += 0.2
        score += 0.4 * min(np.log1p(d['n_interactions']) / np.log1p(90), 1.0)
    if d['n_approved'] > 0:
        score += 0.2
    if d['has_agonist']:
        score += 0.2
    return round(min(score, 1.0), 3)

# Literature score:
# Based on PubMed count for gene + cardiac terms
# Normalize: log-scale
def literature_score(gene):
    count = pubmed_results.get(gene, 0)
    if count == 0:
        return 0.0
    # log normalize: 1 paper = 0.1, 10 = 0.3, 100 = 0.6, 1000 = 0.9
    return round(min(np.log10(count + 1) / np.log10(1500), 1.0), 3)

# Compute for all genes
gene_scores = {}
for gene in all_genes:
    gene_scores[gene] = {
        'druggability_auto': druggability_score(gene),
        'literature_auto': literature_score(gene),
        'pubmed_count': pubmed_results.get(gene, 0),
        'dgidb_interactions': dgidb_results.get(gene, {}).get('n_interactions', 0),
        'dgidb_approved': dgidb_results.get(gene, {}).get('n_approved', 0),
    }

scores_df = pd.DataFrame(gene_scores).T
scores_df.index.name = 'gene'
scores_df = scores_df.sort_values('literature_auto', ascending=False)

print("\nTop 20 genes by automated scores:")
print(scores_df.head(20).to_string())

# ============================================================================
# 5. Re-run prioritization with automated scores
# ============================================================================
print("\n" + "=" * 80)
print("Step 5: Re-running Phase 6 prioritization with automated scores")
print("=" * 80)

# Use mouse refined pairs as base
pairs = mouse_refined[['receptor', 'ligand', 'r_logfc', 'l_logfc',
                        'n_db', 'starvation', 'composite', 'pathway']].copy()
pairs = pairs.rename(columns={'composite': 'mismatch_composite'})

# Map to human gene names for cross-species lookup
pairs['h_receptor'] = pairs['receptor'].str.upper()
pairs['h_ligand'] = pairs['ligand'].str.upper()

# Conservation (from cross-species data)
cross_all = pd.read_csv(MISMATCH_DIR / "cross_species_lr_mismatch.csv")
conservation_map = {}
for _, row in cross_all.iterrows():
    conservation_map[(row['receptor'], row['ligand'])] = row['conservation']

conservation_scores = {
    'CONSERVED_strict': 1.0, 'CONSERVED_relaxed': 0.7,
    'CONSERVED_trend': 0.5, 'mouse_only': 0.2,
    'human_only': 0.1, 'neither': 0.0, 'not_checked': 0.1,
}

pairs['conservation'] = pairs.apply(
    lambda r: conservation_map.get((r['h_receptor'], r['h_ligand']), 'not_checked'), axis=1)
pairs['conservation_score'] = pairs['conservation'].map(conservation_scores)

# Automated druggability: average of receptor + ligand scores
# (for "deliver ligand" strategy, ligand druggability matters more)
pairs['drug_receptor'] = pairs['h_receptor'].map(
    lambda g: gene_scores.get(g, {}).get('druggability_auto', 0.1))
pairs['drug_ligand'] = pairs['h_ligand'].map(
    lambda g: gene_scores.get(g, {}).get('druggability_auto', 0.1))
pairs['druggability_auto'] = 0.3 * pairs['drug_receptor'] + 0.7 * pairs['drug_ligand']

# Automated literature: average of receptor + ligand
pairs['lit_receptor'] = pairs['h_receptor'].map(
    lambda g: gene_scores.get(g, {}).get('literature_auto', 0.0))
pairs['lit_ligand'] = pairs['h_ligand'].map(
    lambda g: gene_scores.get(g, {}).get('literature_auto', 0.0))
pairs['literature_auto'] = (pairs['lit_receptor'] + pairs['lit_ligand']) / 2

# Normalize mismatch
max_m = pairs['mismatch_composite'].max()
pairs['mismatch_norm'] = pairs['mismatch_composite'] / max_m if max_m > 0 else 0

# Final score
W_MISMATCH = 0.30
W_CONSERVATION = 0.30
W_DRUG = 0.20
W_LIT = 0.20

pairs['priority_score'] = (
    W_MISMATCH * pairs['mismatch_norm'] +
    W_CONSERVATION * pairs['conservation_score'] +
    W_DRUG * pairs['druggability_auto'] +
    W_LIT * pairs['literature_auto']
)

pairs = pairs.sort_values('priority_score', ascending=False).reset_index(drop=True)
pairs.index += 1
pairs.index.name = 'rank'

# ============================================================================
# 6. Display results
# ============================================================================
print("\nFINAL RANKING (automated scoring)")
print(f"Weights: Mismatch {W_MISMATCH:.0%} + Conservation {W_CONSERVATION:.0%} "
      f"+ Druggability {W_DRUG:.0%} + Literature {W_LIT:.0%}")

cols = ['receptor', 'ligand', 'priority_score', 'mismatch_norm',
        'conservation', 'conservation_score', 'druggability_auto',
        'literature_auto', 'n_db', 'pathway']
print("\nTop 30:")
print(pairs[cols].head(30).to_string())

# Score breakdown top 15
print(f"\n{'Rank':<5} {'Pair':<22} {'TOTAL':>6} = "
      f"{'Mismatch':>8}  + {'Conserv':>7}  + {'Drug':>6}  + {'Lit':>6}  {'Pathway'}")
print("-" * 95)
for i in range(min(15, len(pairs))):
    r = pairs.iloc[i]
    pair = f"{r['receptor']}/{r['ligand']}"
    print(f" {i+1:<4} {pair:<22} {r['priority_score']:.3f} = "
          f"{W_MISMATCH*r['mismatch_norm']:>8.3f}  + {W_CONSERVATION*r['conservation_score']:>7.3f}  + "
          f"{W_DRUG*r['druggability_auto']:>6.3f}  + {W_LIT*r['literature_auto']:>6.3f}  {r.get('pathway','')}")

# Key pairs
print("\nKey pairs:")
for rec, lig in [('Fgfr2','Fgf10'),('Fgfr2','Fgf7'),('Acvr1','Bmp6'),
                 ('Bmpr2','Bmp4'),('Tyro3','Gas6'),('Unc5b','Ntn1')]:
    m = pairs[(pairs['receptor']==rec) & (pairs['ligand']==lig)]
    if len(m) > 0:
        idx = m.index[0]
        row = m.iloc[0]
        print(f"  {rec}/{lig:<10} rank={idx:>3}/{len(pairs)}  score={row['priority_score']:.3f}  "
              f"drug_auto={row['druggability_auto']:.3f}  lit_auto={row['literature_auto']:.3f}")

# ============================================================================
# 7. Save
# ============================================================================
print("\n" + "=" * 80)
print("Saving")
print("=" * 80)

pairs.to_csv(MISMATCH_DIR / "therapeutic_targets_automated.csv")
scores_df.to_csv(MISMATCH_DIR / "gene_scores_automated.csv")
print(f"Saved: therapeutic_targets_automated.csv ({len(pairs)} pairs)")
print(f"Saved: gene_scores_automated.csv ({len(scores_df)} genes)")
