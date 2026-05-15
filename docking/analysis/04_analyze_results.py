"""
Script 04: Correlate AutoDock Vina docking scores with PDSP pKi values.

Merges docking scores from script 03 with measured 5-HT2A pKi values
from the PDSP Ki database (shulgin-selectivity analysis).

Key questions:
  1. Does docking score correlate with pKi overall?
  2. Do tryptamines and phenethylamines cluster separately?
     (scaffold confounding in 3D — analogous to QSAR finding)
  3. Which compounds are well-predicted vs outliers?

pKi values are loaded from the selectivity profiles computed in
shulgin-selectivity/analysis/selectivity_profiles.csv.

Outputs:
  analysis/04_correlation_report.txt
  analysis/docking_vs_pki.png          - scatter plot coloured by scaffold
  analysis/docking_score_distribution.png - score distribution by scaffold class
  results/docking_pki_merged.csv       - merged data for further analysis

Usage:
  cd psychedelic-selectivity/docking/
  python analysis/04_analyze_results.py
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
import os

log = open("analysis/04_correlation_report.txt", "w")
def p(*args, **kwargs):
    print(*args, **kwargs, file=log)

# ── 1. Load docking scores ────────────────────────────────────────────────────

scores = pd.read_csv("results/docking_scores.csv")
scores = scores[scores["status"] == "OK"][["compound", "best_score"]].copy()
print(f"Docking scores loaded: {len(scores)} compounds")

# ── 2. Load PDSP pKi values ───────────────────────────────────────────────────

PKI_CSV = "../shulgin-selectivity/analysis/selectivity_profiles.csv"

try:
    sel = pd.read_csv(PKI_CSV)
    # Column is named 'compound' and '5-HT2A' (the pKi column)
    # selectivity_profiles.csv has receptor columns as pKi values
    pki = sel[["compound", "5-HT2A"]].rename(columns={"5-HT2A": "pki_2A"})
    pki = pki[pki["pki_2A"].notna()]
    print(f"PDSP pKi values loaded: {len(pki)} compounds")
except FileNotFoundError:
    print(f"ERROR: {PKI_CSV} not found.")
    print("Run shulgin-selectivity analysis first.")
    raise

# ── 3. Merge ──────────────────────────────────────────────────────────────────

df = scores.merge(pki, on="compound", how="inner")
print(f"Merged: {len(df)} compounds with both docking score and pKi")

# Assign scaffold class
def scaffold_class(name):
    tryptamines = {"DMT", "5-MeO-DMT", "bufotenine", "psilocin", "psilocybin",
                   "5-Me-DMT", "6-F-DMT", "DPT", "DiPT", "5-MeO-DiPT",
                   "5-MeO-MiPT", "5-MeO-T", "AMT", "AMT (+)", "AMT (-)",
                   "DALT", "5-MeO-DALT", "5-F-DALT", "5-Br-DALT",
                   "4-HO-DALT", "4-AcO-DALT", "2-Ph-DALT", "7-Et-DALT",
                   "5-MeO-2-Me-DALT", "5-MeO-2-F-DALT",
                   "tryptamine", "serotonin", "melatonin", "ibogaine"}
    phenethylamines = {"mescaline", "MDMA", "MDA", "MDA R(-)", "MDA (R,S)",
                       "DOB", "DOI", "DOM", "DOET",
                       "2C-B", "2C-E", "2C-T-2", "2C-H"}
    ergolines = {"LSD", "LSD (+)", "2-Bromo-LSD"}
    if name in tryptamines:   return "tryptamine"
    if name in phenethylamines: return "phenethylamine"
    if name in ergolines:     return "ergoline"
    return "other"

df["scaffold_class"] = df["compound"].apply(scaffold_class)

df.to_csv("results/docking_pki_merged.csv", index=False)

# ── 4. Correlation analysis ───────────────────────────────────────────────────

p("Docking Score vs PDSP pKi Correlation Analysis")
p("=" * 65)
p(f"Compounds with both docking score and PDSP pKi: {len(df)}")
p(f"\nDocking score range: {df['best_score'].min():.3f} to "
  f"{df['best_score'].max():.3f} kcal/mol")
p(f"pKi range: {df['pki_2A'].min():.3f} to {df['pki_2A'].max():.3f}")

# Overall correlation (note: more negative score = better binding,
# so expect NEGATIVE correlation with pKi if docking is predictive)
r, pval = stats.pearsonr(df["best_score"], df["pki_2A"])
rho, pval_s = stats.spearmanr(df["best_score"], df["pki_2A"])
p(f"\nOverall correlation (n={len(df)}):")
p(f"  Pearson  r  = {r:.3f}  (p={pval:.3e})")
p(f"  Spearman rho = {rho:.3f}  (p={pval_s:.3e})")
p(f"  Note: negative r = docking score predicts binding "
  f"(more negative score -> higher pKi)")

# Per-scaffold correlation
p(f"\nCorrelation within scaffold classes:")
for cls in ["tryptamine", "phenethylamine"]:
    sub = df[df["scaffold_class"] == cls]
    if len(sub) < 4:
        p(f"  {cls}: n={len(sub)} — too few for correlation")
        continue
    r_c, p_c = stats.pearsonr(sub["best_score"], sub["pki_2A"])
    p(f"  {cls:<16} n={len(sub):>2}  r={r_c:.3f}  p={p_c:.3e}")

# Score distribution by scaffold
p(f"\nMean docking score by scaffold class:")
for cls, grp in df.groupby("scaffold_class"):
    p(f"  {cls:<16} n={len(grp):>2}  "
      f"mean={grp['best_score'].mean():.3f}  "
      f"std={grp['best_score'].std():.3f}  "
      f"range=[{grp['best_score'].min():.3f}, {grp['best_score'].max():.3f}]")

# Notable outliers
p(f"\nFull compound table (sorted by docking score):")
p(f"{'Compound':<28} {'Score':>8} {'pKi_2A':>8} {'Class':<16}")
p("-" * 65)
for _, row in df.sort_values("best_score").iterrows():
    p(f"  {row['compound']:<26} {row['best_score']:>8.3f} "
      f"{row['pki_2A']:>8.3f} {row['scaffold_class']:<16}")

# ── 5. Plots ──────────────────────────────────────────────────────────────────

colors = {"tryptamine": "#4C72B0", "phenethylamine": "#DD8452",
          "ergoline": "#55A868", "other": "#888888"}
markers = {"tryptamine": "o", "phenethylamine": "s",
           "ergoline": "^", "other": "D"}

# Plot 1: docking score vs pKi scatter
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

ax = axes[0]
for cls in df["scaffold_class"].unique():
    sub = df[df["scaffold_class"] == cls]
    ax.scatter(sub["best_score"], sub["pki_2A"],
               color=colors[cls], marker=markers[cls],
               s=60, alpha=0.8, label=cls, zorder=3)
    for _, row in sub.iterrows():
        ax.annotate(row["compound"][:12],
                    (row["best_score"], row["pki_2A"]),
                    fontsize=5.5, alpha=0.7,
                    xytext=(3, 3), textcoords="offset points")

# Overall regression line
x_fit = np.linspace(df["best_score"].min(), df["best_score"].max(), 100)
slope, intercept, r_val, p_val, _ = stats.linregress(
    df["best_score"], df["pki_2A"])
ax.plot(x_fit, slope * x_fit + intercept,
        "k--", lw=1, alpha=0.5,
        label=f"Overall fit (r={r_val:.2f}, p={p_val:.3f})")

ax.set_xlabel("AutoDock Vina score (kcal/mol)\n"
              "(more negative = better predicted binding)")
ax.set_ylabel("PDSP pKi at 5-HT2A\n(higher = stronger measured binding)")
ax.set_title("Docking score vs measured binding affinity\n"
             "5-HT2A (6WHA active state), Shulgin compounds")
ax.legend(fontsize=7)

# Per-scaffold regression lines
for cls in ["tryptamine", "phenethylamine"]:
    sub = df[df["scaffold_class"] == cls]
    if len(sub) < 4: continue
    sl, ic, rv, pv, _ = stats.linregress(sub["best_score"], sub["pki_2A"])
    x_s = np.linspace(sub["best_score"].min(), sub["best_score"].max(), 50)
    ax.plot(x_s, sl * x_s + ic, color=colors[cls],
            lw=1.5, linestyle=":", alpha=0.7,
            label=f"{cls} (r={rv:.2f})")

ax.legend(fontsize=7)

# Plot 2: score distribution by scaffold class (box + strip)
ax2 = axes[1]
classes_ordered = ["tryptamine", "phenethylamine"]
data_by_class   = [df[df["scaffold_class"]==c]["best_score"].values
                   for c in classes_ordered]
bp = ax2.boxplot(data_by_class, labels=classes_ordered,
                 patch_artist=True, widths=0.4)
for patch, cls in zip(bp["boxes"], classes_ordered):
    patch.set_facecolor(colors[cls])
    patch.set_alpha(0.6)

# Overlay individual points
for xi, (cls, data) in enumerate(zip(classes_ordered, data_by_class),
                                  start=1):
    jitter = np.random.default_rng(42).uniform(-0.1, 0.1, len(data))
    ax2.scatter(xi + jitter, data, color=colors[cls],
                s=40, alpha=0.7, zorder=3)

# Annotate means
for xi, data in enumerate(data_by_class, start=1):
    ax2.text(xi, np.mean(data) + 0.05, f"μ={np.mean(data):.2f}",
             ha="center", fontsize=8, color="black")

ax2.set_ylabel("AutoDock Vina score (kcal/mol)")
ax2.set_title("Score distribution by scaffold class\n"
              "(separation quantifies scaffold bias in docking)")

# T-test between classes
t, tp = stats.ttest_ind(data_by_class[0], data_by_class[1])
ax2.text(0.5, 0.02, f"t-test: t={t:.2f}, p={tp:.4f}",
         transform=ax2.transAxes, ha="center", fontsize=8,
         bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

plt.suptitle("AutoDock Vina at 5-HT2A (6WHA): docking scores for Shulgin compounds",
             fontsize=10)
plt.tight_layout()
plt.savefig("analysis/docking_vs_pki.png", dpi=150)
p("\nSaved analysis/docking_vs_pki.png")
plt.close()

log.close()
print("Done. Report: analysis/04_correlation_report.txt")
print(f"Overall Pearson r = {r:.3f} (p={pval:.3e})")
print(f"Tryptamines vs phenethylamines: t-test p = {tp:.4f}")
