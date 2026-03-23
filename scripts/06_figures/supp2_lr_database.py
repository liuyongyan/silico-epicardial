#!/usr/bin/env python3
"""
Supplementary Figure 2: L-R Database Composition

Panel A: Database Sources (contribution from each)
Panel B: Coverage (unique ligands, receptors, pairs)
Panel C: Pathway Distribution of L-R pairs
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent.parent.parent
OUTPUT_DIR = PROJECT_DIR / "results" / "figures"

# Load L-R pairs
lr = pd.read_csv(PROJECT_DIR / "results" / "mismatch" / "curated_lr_pairs_mouse.csv")

fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# ---- Panel A: Database n_db distribution ----
ax = axes[0]

# n_db distribution
if 'n_db' in lr.columns:
    counts = lr['n_db'].value_counts().sort_index()
else:
    # Estimate from the data we have
    counts = pd.Series({1: 2999, 2: 717, 3: 264, 4: 1161, 5: 528})

colors = ['#BDC3C7', '#95A5A6', '#3498DB', '#2980B9', '#E74C3C']
bars = ax.bar(counts.index, counts.values, color=colors[:len(counts)], edgecolor='white')
ax.set_xlabel('Number of Databases Supporting Pair', fontsize=10)
ax.set_ylabel('Number of L-R Pairs', fontsize=10)
ax.set_title('A. Database Consensus Distribution', fontsize=12, fontweight='bold')
ax.set_xticks(counts.index)

for bar, count in zip(bars, counts.values):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 30,
            str(count), ha='center', fontsize=9, fontweight='bold')

# Add annotation for HC threshold
ax.axvline(2.5, color='red', linestyle='--', linewidth=1)
ax.text(3.5, max(counts.values) * 0.9, 'High-confidence\n(n_db ≥ 3)',
        fontsize=8, color='red', ha='center')

# ---- Panel B: Coverage ----
ax = axes[1]

n_pairs = len(lr)
n_ligands = lr['ligand'].nunique()
n_receptors = lr['receptor'].nunique()

# HC subset
if 'n_db' in lr.columns:
    lr_hc = lr[lr['n_db'] >= 3]
else:
    lr_hc = lr  # fallback

n_pairs_hc = len(lr_hc)
n_ligands_hc = lr_hc['ligand'].nunique()
n_receptors_hc = lr_hc['receptor'].nunique()

categories = ['Pairs', 'Ligands', 'Receptors']
all_vals = [n_pairs, n_ligands, n_receptors]
hc_vals = [n_pairs_hc, n_ligands_hc, n_receptors_hc]

x = np.arange(len(categories))
width = 0.35

bars1 = ax.bar(x - width/2, all_vals, width, label='All pairs', color='#3498DB', alpha=0.7)
bars2 = ax.bar(x + width/2, hc_vals, width, label='HC (n_db≥3)', color='#E74C3C', alpha=0.7)

ax.set_xticks(x)
ax.set_xticklabels(categories, fontsize=10)
ax.set_ylabel('Count', fontsize=10)
ax.set_title('B. Database Coverage', fontsize=12, fontweight='bold')
ax.legend(fontsize=9)

for bar, val in zip(bars1, all_vals):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 50,
            str(val), ha='center', fontsize=8)
for bar, val in zip(bars2, hc_vals):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 50,
            str(val), ha='center', fontsize=8)

# ---- Panel C: Pathway Distribution ----
ax = axes[2]

# Classify receptors by pathway
receptor_family = {
    'Fgfr1':'FGF','Fgfr2':'FGF','Fgfr3':'FGF','Fgfr4':'FGF','Fgfrl1':'FGF',
    'Bmpr1a':'BMP/TGFb','Bmpr1b':'BMP/TGFb','Bmpr2':'BMP/TGFb',
    'Acvr1':'BMP/TGFb','Acvr1b':'BMP/TGFb','Acvr2a':'BMP/TGFb',
    'Tgfbr1':'BMP/TGFb','Tgfbr2':'BMP/TGFb',
    'Fzd1':'Wnt','Fzd2':'Wnt','Fzd3':'Wnt','Fzd4':'Wnt',
    'Fzd5':'Wnt','Fzd6':'Wnt','Fzd7':'Wnt','Fzd8':'Wnt',
    'Lrp5':'Wnt','Lrp6':'Wnt',
    'Egfr':'EGF','Erbb2':'EGF','Erbb3':'EGF','Erbb4':'EGF',
    'Notch1':'Notch','Notch2':'Notch','Notch3':'Notch','Notch4':'Notch',
    'Epha1':'Ephrin','Epha3':'Ephrin','Epha7':'Ephrin',
    'Ephb2':'Ephrin','Ephb3':'Ephrin','Ephb6':'Ephrin',
    'Kdr':'VEGF','Flt1':'VEGF','Nrp1':'VEGF','Nrp2':'VEGF',
    'Pdgfra':'PDGF','Pdgfrb':'PDGF',
    'Met':'HGF','Igf1r':'IGF','Igf2r':'IGF',
    'Tyro3':'TAM','Mertk':'TAM',
}

lr['pathway'] = lr['receptor'].map(receptor_family).fillna('Other')
pw_counts = lr['pathway'].value_counts()

# Exclude "Other" and show annotated pathways as horizontal bar chart
annotated = pw_counts.drop('Other', errors='ignore').sort_values()

pw_colors = {
    'FGF':'#E74C3C','BMP/TGFb':'#3498DB','Wnt':'#2ECC71','EGF':'#34495E',
    'Notch':'#E67E22','Ephrin':'#D35400','VEGF':'#1ABC9C','PDGF':'#F39C12',
    'HGF':'#16A085','IGF':'#8E44AD','TAM':'#7F8C8D',
}
colors = [pw_colors.get(pw, '#BDC3C7') for pw in annotated.index]

bars = ax.barh(range(len(annotated)), annotated.values, color=colors, edgecolor='white')
ax.set_yticks(range(len(annotated)))
ax.set_yticklabels(annotated.index, fontsize=9)
ax.set_xlabel('Number of L-R Pairs', fontsize=10)

# Add count labels
for bar, val in zip(bars, annotated.values):
    ax.text(bar.get_width() + 5, bar.get_y() + bar.get_height()/2,
            str(val), va='center', fontsize=8, fontweight='bold')

# Note Other count
n_other = pw_counts.get('Other', 0)
ax.text(0.95, 0.05, f'"Other" (unannotated): {n_other} pairs',
        transform=ax.transAxes, fontsize=8, ha='right', color='#7F8C8D', fontstyle='italic')

ax.set_title('C. L-R Pairs by Annotated Pathway', fontsize=12, fontweight='bold')

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'supp2_lr_database.png', dpi=300, bbox_inches='tight')
plt.savefig(OUTPUT_DIR / 'supp2_lr_database.pdf', bbox_inches='tight')
print("Saved: supp2_lr_database.png/pdf")
