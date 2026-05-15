"""
Script 05: Summary figure — docking scaffold bias at 5-HT2A.

Produces a three-panel figure connecting the docking results to the
broader scaffold confounding narrative:

Panel 1: Docking score vs measured pKi, coloured by scaffold class.
         Shows the overall positive (wrong-direction) correlation and
         the complete separation between tryptamines and phenethylamines.

Panel 2: Docking score distribution by scaffold class (box + strip).
         Shows tryptamines systematically score better despite comparable
         or lower measured binding affinity.

Panel 3: Within-class docking score vs pKi.
         Shows phenethylamines have r=+0.80 in the wrong direction —
         the strongest binders (DOB, DOI, 2C-B) dock worst.

Outputs:
  analysis/docking_summary.png
  analysis/05_summary_report.txt

Usage:
  cd psychedelic-selectivity/docking/
  python analysis/05_summary_figure.py
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from scipy import stats
import os

log = open("analysis/05_summary_report.txt", "w")
def p(*args, **kwargs):
    print(*args, **kwargs, file=log)

# ── Load merged data ──────────────────────────────────────────────────────────

df = pd.read_csv("results/docking_pki_merged.csv")
df = df[df["best_score"].notna() & df["pki_2A"].notna()]

trp  = df[df["scaffold_class"] == "tryptamine"]
phe  = df[df["scaffold_class"] == "phenethylamine"]

colors  = {"tryptamine": "#4C72B0", "phenethylamine": "#DD8452",
           "ergoline": "#55A868", "other": "#888888"}
markers = {"tryptamine": "o", "phenethylamine": "s",
           "ergoline": "^", "other": "D"}

# ── Compute statistics ────────────────────────────────────────────────────────

r_all,   p_all   = stats.pearsonr(df["best_score"],  df["pki_2A"])
r_trp,   p_trp   = stats.pearsonr(trp["best_score"], trp["pki_2A"])
r_phe,   p_phe   = stats.pearsonr(phe["best_score"], phe["pki_2A"])
rho_all, p_rho   = stats.spearmanr(df["best_score"], df["pki_2A"])
t_stat,  t_pval  = stats.ttest_ind(trp["best_score"], phe["best_score"])

p("Docking Summary Report")
p("=" * 65)
p(f"n = {len(df)} compounds ({len(trp)} tryptamines, {len(phe)} phenethylamines)")
p(f"\nOverall:        Pearson r = {r_all:+.3f}  p = {p_all:.4f}")
p(f"                Spearman rho = {rho_all:+.3f}  p = {p_rho:.4f}")
p(f"\nTryptamines:    r = {r_trp:+.3f}  p = {p_trp:.4f}  (n={len(trp)})")
p(f"Phenethylamines:r = {r_phe:+.3f}  p = {p_phe:.4f}  (n={len(phe)})")
p(f"\nScaffold score gap: tryptamine mean = {trp['best_score'].mean():.3f}, "
  f"phenethylamine mean = {phe['best_score'].mean():.3f}")
p(f"  t-test: t = {t_stat:.2f}, p = {t_pval:.2e}")
p(f"\nInterpretation:")
p(f"  Overall r = {r_all:+.3f} is positive — WRONG direction.")
p(f"  A predictive model should have r < 0 (better docking = higher pKi).")
p(f"  The positive correlation is driven by scaffold separation:")
p(f"    tryptamines dock better (mean {trp['best_score'].mean():.2f})")
p(f"    phenethylamines dock worse (mean {phe['best_score'].mean():.2f})")
p(f"  but phenethylamines have comparable or higher measured pKi.")
p(f"  DOB (pKi=8.27) and DOI (pKi=8.28) are among the strongest binders")
p(f"  yet dock at only -5.57 and -5.66 kcal/mol.")
p(f"\n  Within phenethylamines: r = {r_phe:+.3f} — strong positive correlation,")
p(f"  again wrong direction. The best binders (2C-B, DOB, DOI) dock worst.")
p(f"\n  Cause: 6WHA was co-crystallised with 25CN-NBOH, a tryptamine-like")
p(f"  compound. The orthosteric pocket geometry is optimised for tryptamine")
p(f"  scaffolds. Phenethylamines bind through a different mode not captured")
p(f"  by this structure.")
p(f"\n  This is the 3D analogue of the QSAR scaffold confounding finding:")
p(f"  structure-based virtual screening at 6WHA will systematically favour")
p(f"  tryptamines over phenethylamines regardless of true affinity.")

# ── Three-panel summary figure ────────────────────────────────────────────────

fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# ── Panel 1: Overall scatter ──────────────────────────────────────────────────

ax1 = axes[0]
for cls in ["tryptamine", "phenethylamine"]:
    sub = df[df["scaffold_class"] == cls]
    ax1.scatter(sub["best_score"], sub["pki_2A"],
                color=colors[cls], marker=markers[cls],
                s=55, alpha=0.85, label=cls, zorder=3)
    for _, row in sub.iterrows():
        ax1.annotate(row["compound"][:10],
                     (row["best_score"], row["pki_2A"]),
                     fontsize=5, alpha=0.65,
                     xytext=(3, 2), textcoords="offset points")

# Overall regression
x_fit  = np.linspace(df["best_score"].min() - 0.1,
                      df["best_score"].max() + 0.1, 100)
sl, ic = np.polyfit(df["best_score"], df["pki_2A"], 1)
ax1.plot(x_fit, sl * x_fit + ic, "k--", lw=1.2, alpha=0.6,
         label=f"All: r={r_all:+.2f} (p={p_all:.3f})")

ax1.set_xlabel("Vina score (kcal/mol)\n(more negative = better predicted binding)",
               fontsize=8)
ax1.set_ylabel("PDSP pKi at 5-HT2A\n(higher = stronger measured binding)",
               fontsize=8)
ax1.set_title("Docking score vs measured pKi\n"
              "positive r = wrong direction", fontsize=9)
ax1.legend(fontsize=7)

# Annotate notable outliers
for _, row in df.iterrows():
    if row["compound"] in ("DOB", "DOI", "2C-B", "2-Ph-DALT", "psilocybin"):
        ax1.annotate(row["compound"],
                     (row["best_score"], row["pki_2A"]),
                     fontsize=6.5, color="black", fontweight="bold",
                     xytext=(4, -8), textcoords="offset points")

# ── Panel 2: Score distribution by class ─────────────────────────────────────

ax2 = axes[1]
classes_ordered = ["tryptamine", "phenethylamine"]
data_by_class   = [df[df["scaffold_class"]==c]["best_score"].values
                   for c in classes_ordered]

bp = ax2.boxplot(data_by_class, tick_labels=classes_ordered,
                 patch_artist=True, widths=0.45)
for patch, cls in zip(bp["boxes"], classes_ordered):
    patch.set_facecolor(colors[cls])
    patch.set_alpha(0.5)

rng = np.random.default_rng(42)
for xi, (cls, data) in enumerate(zip(classes_ordered, data_by_class), 1):
    jitter = rng.uniform(-0.12, 0.12, len(data))
    ax2.scatter(xi + jitter, data, color=colors[cls], s=40, alpha=0.8, zorder=4)

for xi, (cls, data) in enumerate(zip(classes_ordered, data_by_class), 1):
    ax2.text(xi, np.mean(data) + 0.05, f"μ={np.mean(data):.2f}",
             ha="center", fontsize=8.5, fontweight="bold")

ax2.set_ylabel("Vina score (kcal/mol)", fontsize=8)
ax2.set_title(f"Score gap between scaffold classes\n"
              f"(t={t_stat:.1f}, p={t_pval:.1e} — complete separation)",
              fontsize=9)
ax2.text(0.5, 0.04,
         f"Tryptamines dock better by {abs(trp['best_score'].mean() - phe['best_score'].mean()):.2f} kcal/mol\n"
         f"despite similar measured affinity",
         transform=ax2.transAxes, ha="center", fontsize=7.5,
         bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))

# ── Panel 3: Within-class correlation ────────────────────────────────────────

ax3 = axes[2]
for cls, r_c, p_c in [("tryptamine", r_trp, p_trp),
                       ("phenethylamine", r_phe, p_phe)]:
    sub = df[df["scaffold_class"] == cls]
    ax3.scatter(sub["best_score"], sub["pki_2A"],
                color=colors[cls], marker=markers[cls],
                s=55, alpha=0.85, label=f"{cls} (r={r_c:+.2f})", zorder=3)
    sl_c, ic_c = np.polyfit(sub["best_score"], sub["pki_2A"], 1)
    x_c = np.linspace(sub["best_score"].min()-0.05,
                       sub["best_score"].max()+0.05, 50)
    ax3.plot(x_c, sl_c * x_c + ic_c, color=colors[cls],
             lw=1.5, linestyle="--", alpha=0.8)
    for _, row in sub.iterrows():
        if row["compound"] in ("DOB","DOI","2C-B","psilocybin","2-Ph-DALT",
                                "5-MeO-DiPT","DALT"):
            ax3.annotate(row["compound"][:9],
                         (row["best_score"], row["pki_2A"]),
                         fontsize=6, color=colors[cls],
                         xytext=(4, 2), textcoords="offset points")

ax3.set_xlabel("Vina score (kcal/mol)", fontsize=8)
ax3.set_ylabel("PDSP pKi at 5-HT2A", fontsize=8)
ax3.set_title("Within-class correlation\n"
              "phenethylamine r=+0.80: best binders dock worst", fontsize=9)
ax3.legend(fontsize=7)

ax3.annotate("DOB pKi=8.27\nbut score=-5.57",
             xy=(-5.572, 8.269), fontsize=7, color="#DD8452",
             xytext=(-7.0, 7.8),
             arrowprops=dict(arrowstyle="->", color="#DD8452", lw=1.0))

plt.suptitle(
    "AutoDock Vina at 5-HT2A (6WHA): scaffold-dependent docking performance\n"
    "Tryptamines systematically score better than phenethylamines "
    "despite comparable measured affinity — 3D scaffold confounding",
    fontsize=9.5)
plt.tight_layout()
plt.savefig("analysis/docking_summary.png", dpi=150, bbox_inches="tight")
p("\nSaved analysis/docking_summary.png")

log.close()
print("Done.")
print(f"Overall r = {r_all:+.3f}  (tryptamine r={r_trp:+.3f}, "
      f"phenethylamine r={r_phe:+.3f})")
print("Figure: analysis/docking_summary.png")
