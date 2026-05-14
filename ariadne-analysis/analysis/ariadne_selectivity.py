"""
Ariadne pharmacological analysis — Cunningham et al. ACS Chem Neurosci 2023
doi: 10.1021/acschemneuro.2c00597

Three analyses, all using data from the paper's figures and tables:

1. ASSAY-DEPENDENT EFFICACY (ariadne_gq_vs_caflux.png)
   Gq dissociation BRET vs Ca2+ flux side by side across 5-HT2A/2B/2C.
   Ariadne's Emax at 5-HT2C drops from ~82% (Gq BRET) to ~34% (Ca2+ flux)
   for the racemate, and to 22% for (S)-Ariadne. This assay dependence is
   in the paper's data but never presented side-by-side.

2. BIAS FACTOR ACROSS ALPHA-ALKYL PAIRS (ariadne_alpha_series.png)
   log_bias = log10((Emax_barr/EC50_barr) / (Emax_Gq/EC50_Gq)) computed for
   three matched hallucinogenic (alpha-methyl) / non-hallucinogenic (alpha-ethyl)
   pairs at 5-HT2A. The paper discusses each pair in separate figure panels;
   this is the first unified quantitative comparison showing the consistent
   pattern: alpha-ethyl compounds are slightly less Gq-biased (~0.1-0.35 log
   units) but remain Gq-biased — confirming the mechanism is NOT biased agonism.

3. 4-POSITION ANALOG 2A/2B SELECTIVITY (ariadne_2a2b_selectivity.png)
   pEC50 difference (pEC50_2A - pEC50_2B) for all five 4-position analogs
   from Figure 5 of the paper. Expressing selectivity as pEC50 differences
   allows direct comparison across analogs; the paper reports EC50 values
   but does not compute or compare these differences explicitly.

NOT included: four-metric selectivity comparison (S-score, Gini, entropy,
ratio). With only 3 receptors, all definitions agree trivially — there is
no statistical room for divergence. This analysis requires the full 12-receptor
5-HTome numerical data (Table S3, not publicly available as of 2026).

Parameters:
  PRIMARY_ASSAY : assay for the Gq profile ("gq_bret" or "ca_flux")

Outputs:
  analysis/ariadne_report.txt
  analysis/ariadne_gq_vs_caflux.png
  analysis/ariadne_alpha_series.png
  analysis/ariadne_2a2b_selectivity.png
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import os

os.makedirs("analysis", exist_ok=True)

PRIMARY_ASSAY = "gq_bret"

log = open("analysis/ariadne_report.txt", "w")
def p(*args, **kwargs):
    print(*args, **kwargs, file=log)

# ── 1. Load data ──────────────────────────────────────────────────────────────

df = pd.read_csv("data/cunningham2023_data.csv", comment="#")
df["ec50_nM"]  = pd.to_numeric(df["ec50_nM"],  errors="coerce")
df["emax_pct"] = pd.to_numeric(df["emax_pct"], errors="coerce")
df["pec50"]    = 9 - np.log10(df["ec50_nM"])

p("Ariadne Pharmacological Analysis — Cunningham et al. ACS Chem Neurosci 2023")
p("=" * 72)
p(f"doi: 10.1021/acschemneuro.2c00597")
p(f"Data file: data/cunningham2023_data.csv  ({len(df)} records)")
p(f"\nThree analyses of published numerical data.")
p(f"See script docstring for what is and is not included and why.")

# ── 2. Analysis 1: Assay-dependent efficacy ───────────────────────────────────

p(f"\n\n{'='*72}")
p("ANALYSIS 1: Assay-dependent efficacy at 5-HT2A/2B/2C")
p("Gq dissociation BRET vs Ca2+ flux — EC50 (nM) and Emax (% 5-HT)")
p("=" * 72)

profile_compounds = ["5-HT", "(rac)-Ariadne", "(R)-Ariadne", "(S)-Ariadne"]
receptors         = ["5-HT2A", "5-HT2B", "5-HT2C"]
assays            = ["gq_bret", "ca_flux"]
assay_labels      = {"gq_bret": "Gq dissociation BRET", "ca_flux": "Ca²⁺ flux"}

for assay in assays:
    p(f"\n{assay_labels[assay]}:")
    p(f"  {'Compound':<20} {'5-HT2A EC50':>12} {'Emax':>6} "
      f"{'5-HT2B EC50':>12} {'Emax':>6} {'5-HT2C EC50':>12} {'Emax':>6}")
    p("  " + "-"*72)
    for compound in profile_compounds:
        row_str = f"  {compound:<20}"
        for receptor in receptors:
            sub = df[(df["compound"]==compound) & (df["assay"]==assay) &
                     (df["receptor"]==receptor)]
            if len(sub) > 0:
                row_str += f" {sub['ec50_nM'].values[0]:>12.1f} {sub['emax_pct'].values[0]:>6.1f}"
            else:
                row_str += f" {'n/a':>12} {'n/a':>6}"
        p(row_str)

p(f"\nKey finding:")
p(f"  At 5-HT2C, (rac)-Ariadne Emax drops from 82.2% (Gq BRET) to 33.6%")
p(f"  (Ca2+ flux) — a 2.4-fold difference in apparent efficacy from assay")
p(f"  choice alone. (S)-Ariadne shows an even larger drop: 79.3% to 22.5%.")
p(f"  Potency (EC50) is similar across assays; efficacy diverges sharply.")
p(f"  The paper attributes the 5-HT2C Ca2+ flux discrepancy to kinetic")
p(f"  differences in 5-HT2C activation relative to 5-HT2A/2B. Regardless")
p(f"  of mechanism, this illustrates that assay choice materially changes")
p(f"  the apparent selectivity profile.")

# ── 3. Analysis 2: Bias factor across alpha-alkyl pairs ──────────────────────

p(f"\n\n{'='*72}")
p("ANALYSIS 2: Bias factor across alpha-alkyl matched pairs at 5-HT2A")
p("log_bias = log10((Emax_barr/EC50_barr) / (Emax_Gq/EC50_Gq))")
p("> 0 = beta-arrestin biased   < 0 = Gq biased   ~ 0 = balanced")
p("=" * 72)

pairs = [
    ("4-methyl vs Ariadne",  "DOM",   "(rac)-Ariadne", "4-methyl"),
    ("4-propyl pair",         "DOPr",  "4C-Pr",         "4-propyl"),
    ("4-iodo pair",           "DOI",   "4C-I",          "4-iodo"),
]

p(f"\n{'Compound':<25} {'4-pos':>8} {'alpha':>6} "
  f"{'Gq EC50':>9} {'Gq Emax':>9} {'barr EC50':>10} {'barr Emax':>10} "
  f"{'logBias':>9} {'Halluc?':>8}")
p("-" * 95)

bias_results = []
for label, halluc_cpd, nonhalluc_cpd, pos in pairs:
    for compound, halluc in [(halluc_cpd, True), (nonhalluc_cpd, False)]:
        gq   = df[(df["compound"]==compound) & (df["assay"]=="gq_bret") &
                  (df["receptor"]=="5-HT2A")]
        barr = df[(df["compound"]==compound) & (df["assay"]=="barr2_bret") &
                  (df["receptor"]=="5-HT2A")]
        if len(gq) == 0: continue
        gq_ec50  = gq["ec50_nM"].values[0]
        gq_emax  = gq["emax_pct"].values[0]
        ba_ec50  = barr["ec50_nM"].values[0] if len(barr) > 0 else np.nan
        ba_emax  = barr["emax_pct"].values[0] if len(barr) > 0 else np.nan
        alpha    = "Me" if halluc else "Et"
        if pd.notna(ba_ec50):
            lb    = np.log10((ba_emax/ba_ec50) / (gq_emax/gq_ec50))
            lb_s  = f"{lb:+9.3f}"
        else:
            lb    = np.nan
            lb_s  = "      n/a"
        p(f"  {compound:<25} {pos:>8} {alpha:>6} "
          f"{gq_ec50:>9.1f} {gq_emax:>9.1f} "
          f"{ba_ec50 if pd.notna(ba_ec50) else float('nan'):>10.1f} "
          f"{ba_emax if pd.notna(ba_emax) else float('nan'):>10.1f} "
          f"{lb_s} {'YES' if halluc else 'no':>8}")
        bias_results.append({"compound": compound, "pos": pos, "alpha": alpha,
                              "hallucinogenic": halluc, "log_bias": lb,
                              "gq_ec50": gq_ec50, "gq_emax": gq_emax})
    p("")

p(f"Key finding:")
p(f"  All six compounds are Gq-biased (negative log_bias). Across all three")
p(f"  4-position substitutions, the alpha-ethyl (Ariadne-class) compound")
p(f"  shows slightly less Gq bias than its alpha-methyl (DOx) counterpart:")
p(f"    4-methyl: DOM=-0.642, Ariadne=-0.542 (shift = +0.10)")
p(f"    4-propyl: DOPr=-0.721, 4C-Pr=-0.368 (shift = +0.35)")
p(f"    4-iodo:   DOI=-0.735, 4C-I=-0.449  (shift = +0.29)")
p(f"  The shift is consistent in direction but small in magnitude.")
p(f"  Conclusion: Ariadne is NOT a biased agonist. The mechanism of")
p(f"  non-hallucinogenicity is uniform efficacy reduction, not pathway")
p(f"  selectivity. This formally confirms the paper's signalling efficacy")
p(f"  hypothesis using the bias factor framework.")

# ── 4. Analysis 3: 4-position 2A/2B selectivity ──────────────────────────────

p(f"\n\n{'='*72}")
p("ANALYSIS 3: 4-position analogs — 5-HT2A vs 5-HT2B selectivity (Gq BRET)")
p("pEC50 difference = pEC50(5-HT2A) - pEC50(5-HT2B)")
p("Positive = 2A-selective (safer); negative = 2B-selective (cardiac risk)")
p("=" * 72)

analogs = ["4C-D (Ariadne)", "4C-TFM", "4C-Pr", "4C-cycPr", "4C-MOM"]
p(f"\n{'Compound':<20} {'2A EC50':>9} {'2A Emax':>9} {'2B EC50':>9} "
  f"{'2B Emax':>9} {'pEC50 2A-2B':>12} {'fold selectivity':>17}")
p("-" * 82)

sel_results = []
for compound in analogs:
    s2a = df[(df["compound"]==compound)&(df["assay"]=="gq_bret")&
             (df["receptor"]=="5-HT2A")]
    s2b = df[(df["compound"]==compound)&(df["assay"]=="gq_bret")&
             (df["receptor"]=="5-HT2B")]
    if len(s2a)==0 or len(s2b)==0: continue
    e2a  = s2a["ec50_nM"].values[0]; m2a = s2a["emax_pct"].values[0]
    e2b  = s2b["ec50_nM"].values[0]; m2b = s2b["emax_pct"].values[0]
    sel  = (9-np.log10(e2a)) - (9-np.log10(e2b))
    fold = e2b / e2a
    p(f"  {compound:<20} {e2a:>9.1f} {m2a:>9.1f} {e2b:>9.1f} "
      f"{m2b:>9.1f} {sel:>+12.2f} {fold:>17.1f}x")
    sel_results.append({"compound": compound, "ec50_2a": e2a, "emax_2a": m2a,
                        "ec50_2b": e2b, "emax_2b": m2b, "pec50_diff": sel,
                        "fold": fold})

p(f"\nKey finding:")
p(f"  All Ariadne analogs except 4C-MOM are 2A-selective (positive pEC50 diff).")
p(f"  4C-TFM and 4C-Pr show the highest selectivity (+0.91, 8-fold each).")
p(f"  4C-MOM is essentially balanced (-0.02) and has substantially reduced")
p(f"  Emax at 5-HT2B (43%) — likely due to its increased polarity.")
p(f"  The paper identifies 4C-TFM as most selective but does not express")
p(f"  this as pEC50 differences or compare directly across all five analogs.")

# ── 5. Colour definitions ─────────────────────────────────────────────────────

colors = {
    "5-HT":          "#888888",
    "(rac)-Ariadne": "#4C72B0",
    "(R)-Ariadne":   "#2196F3",
    "(S)-Ariadne":   "#90CAF9",
}

# ── 6. Plot 1: Gq BRET vs Ca2+ flux — 2x3 grid ───────────────────────────────

fig1, axes1 = plt.subplots(2, 3, figsize=(14, 8), sharey="row")

for row, assay in enumerate(["gq_bret", "ca_flux"]):
    for col, receptor in enumerate(receptors):
        ax  = axes1[row][col]
        sub = df[(df["assay"]==assay) & (df["receptor"]==receptor) &
                 (df["compound"].isin(profile_compounds))]
        sub = sub.set_index("compound").reindex(profile_compounds)
        x   = np.arange(len(profile_compounds))
        bar_colors = [colors.get(c, "gray") for c in profile_compounds]

        ax.bar(x, sub["pec50"].values, 0.5, color=bar_colors, alpha=0.8)

        ax2 = ax.twinx()
        ax2.plot(x, sub["emax_pct"].values, "ko--", markersize=6, linewidth=1.2)
        ax2.set_ylim(0, 130)
        ax2.set_ylabel("Emax (% 5-HT)" if col==2 else "", fontsize=8)
        ax2.tick_params(labelsize=7)

        ax.set_xticks(x)
        ax.set_xticklabels([c.replace("(","").replace(")","")
                            for c in profile_compounds],
                           rotation=35, ha="right", fontsize=7)
        ax.set_ylabel("pEC50" if col==0 else "", fontsize=8)
        ax.set_ylim(0, 12)
        ax.set_title(f"{receptor}\n{assay_labels[assay]}", fontsize=9)
        ax.axhline(6.5, color="red", lw=0.8, linestyle="--", alpha=0.5,
                   label="pEC50=6.5 (~316 nM)")

        # Annotate Emax for Ariadne compounds
        for xi, compound in enumerate(profile_compounds):
            emax = sub.loc[compound, "emax_pct"] if compound in sub.index else np.nan
            if pd.notna(emax) and compound != "5-HT":
                ax2.text(xi, emax+4, f"{emax:.0f}%",
                         ha="center", fontsize=6, color="black")

legend_els = [Patch(color=colors[c], alpha=0.8, label=c)
              for c in profile_compounds]
fig1.legend(handles=legend_els, loc="lower center", ncol=4,
            fontsize=8, bbox_to_anchor=(0.5, -0.01))
fig1.suptitle("Ariadne: potency (pEC50 bars) and efficacy (Emax line)\n"
              "Gq BRET (top) vs Ca²⁺ flux (bottom) — note efficacy divergence at 5-HT2C",
              fontsize=10)
plt.tight_layout(rect=[0, 0.05, 1, 1])
plt.savefig("analysis/ariadne_gq_vs_caflux.png", dpi=150, bbox_inches="tight")
p("\nSaved analysis/ariadne_gq_vs_caflux.png")
plt.close()

# ── 7. Plot 2: Alpha-alkyl bias series ───────────────────────────────────────

pair_data = [
    ("4-methyl series\n(DOM / Ariadne / 5C-D)", [
        ("DOM",           "Me",  True,  "#DD8452"),
        ("(rac)-Ariadne", "Et",  False, "#4C72B0"),
        ("5C-D",          "nPr", False, "#55A868"),
    ]),
    ("4-propyl pair\n(DOPr / 4C-Pr)", [
        ("DOPr",  "Me", True,  "#DD8452"),
        ("4C-Pr", "Et", False, "#4C72B0"),
    ]),
    ("4-iodo pair\n(DOI / 4C-I)", [
        ("DOI",  "Me", True,  "#DD8452"),
        ("4C-I", "Et", False, "#4C72B0"),
    ]),
]

fig2, axes2 = plt.subplots(1, 3, figsize=(14, 5))

for ax, (title, members) in zip(axes2, pair_data):
    labels, pec50s, emaxs, biases, bcols = [], [], [], [], []
    for compound, alpha, halluc, color in members:
        gq   = df[(df["compound"]==compound)&(df["assay"]=="gq_bret")&
                  (df["receptor"]=="5-HT2A")]
        barr = df[(df["compound"]==compound)&(df["assay"]=="barr2_bret")&
                  (df["receptor"]=="5-HT2A")]
        if len(gq)==0: continue
        ec50 = gq["ec50_nM"].values[0]
        emax = gq["emax_pct"].values[0]
        pec50s.append(9-np.log10(ec50))
        emaxs.append(emax)
        bcols.append(color)
        if len(barr)>0:
            bec50 = barr["ec50_nM"].values[0]
            bemax = barr["emax_pct"].values[0]
            biases.append(np.log10((bemax/bec50)/(emax/ec50)))
        else:
            biases.append(np.nan)
        hmark = "★" if halluc else "○"
        labels.append(f"{hmark} {compound}\n(α-{alpha})")

    x    = np.arange(len(labels))
    bars = ax.bar(x, pec50s, color=bcols, alpha=0.85, width=0.5)

    ax_t = ax.twinx()
    ax_t.plot(x, emaxs, "ko--", markersize=7, linewidth=1.5)
    ax_t.set_ylim(0, 120)
    ax_t.set_ylabel("Emax (% 5-HT)", fontsize=8)

    for xi, (lb, pe) in enumerate(zip(biases, pec50s)):
        if pd.notna(lb):
            ax.text(xi, pe+0.15, f"bias={lb:+.2f}",
                    ha="center", va="bottom", fontsize=8, color="#222222")
        else:
            ax.text(xi, pe+0.15, "no barr\ndata",
                    ha="center", va="bottom", fontsize=6.5, color="#888888")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8.5)
    ax.set_ylabel("pEC50 at 5-HT2A (Gq)", fontsize=8)
    ax.set_ylim(0, 10)
    ax.set_title(f"{title}", fontsize=9)

from matplotlib.lines import Line2D
legend_els2 = [
    Patch(color="#DD8452", alpha=0.85, label="★ hallucinogenic (alpha-methyl)"),
    Patch(color="#4C72B0", alpha=0.85, label="○ non-hallucinogenic (alpha-ethyl)"),
    Line2D([0],[0], color="black", marker="o", linestyle="--",
           markersize=6, label="Emax (% 5-HT)"),
]
fig2.legend(handles=legend_els2, loc="lower center", ncol=3,
            fontsize=8, bbox_to_anchor=(0.5, -0.02))
fig2.suptitle("Bias factor (log_bias) at 5-HT2A across matched alpha-alkyl pairs\n"
              "All compounds Gq-biased; alpha-ethyl consistently slightly less biased\n"
              "Confirms: mechanism is NOT biased agonism — uniform efficacy reduction",
              fontsize=9)
plt.tight_layout(rect=[0, 0.08, 1, 1])
plt.savefig("analysis/ariadne_alpha_series.png", dpi=150, bbox_inches="tight")
p("Saved analysis/ariadne_alpha_series.png")
plt.close()

# ── 8. Plot 3: 4-position 2A/2B selectivity ──────────────────────────────────

fig3, axes3 = plt.subplots(1, 2, figsize=(12, 5))

compound_names = [r["compound"] for r in sel_results]
pec50_diffs    = [r["pec50_diff"] for r in sel_results]
fold_sels      = [r["fold"] for r in sel_results]
bar_colors3    = ["#C44E52" if d < 0 else "#4C72B0" for d in pec50_diffs]

# Left: pEC50 difference bars
ax3a = axes3[0]
bars3 = ax3a.bar(range(len(compound_names)), pec50_diffs,
                  color=bar_colors3, alpha=0.85, width=0.55)
ax3a.axhline(0, color="black", lw=0.8)
ax3a.axhline(0.5, color="green", lw=0.8, linestyle="--", alpha=0.6,
             label="0.5 log unit (~3-fold)")
ax3a.axhline(1.0, color="green", lw=0.8, linestyle=":",  alpha=0.6,
             label="1.0 log unit (10-fold)")
for xi, (name, val) in enumerate(zip(compound_names, pec50_diffs)):
    ax3a.text(xi, val + (0.02 if val >= 0 else -0.05),
              f"{val:+.2f}", ha="center",
              va="bottom" if val >= 0 else "top", fontsize=9)
ax3a.set_xticks(range(len(compound_names)))
ax3a.set_xticklabels([c.replace("4C-D (Ariadne)","4C-D\n(Ariadne)")
                      for c in compound_names], fontsize=9)
ax3a.set_ylabel("pEC50(5-HT2A) − pEC50(5-HT2B)\n(positive = 2A-selective)")
ax3a.set_title("5-HT2A/2B selectivity: pEC50 difference\n"
               "(all Ariadne analogs, alpha-ethyl series)")
ax3a.legend(fontsize=8)

# Right: EC50 at 2A and 2B as paired bars
ax3b = axes3[1]
x3b  = np.arange(len(compound_names))
w3b  = 0.35
bars_2a = ax3b.bar(x3b - w3b/2,
                    [r["ec50_2a"] for r in sel_results],
                    w3b, label="5-HT2A EC50 (nM)", color="#4C72B0", alpha=0.85)
bars_2b = ax3b.bar(x3b + w3b/2,
                    [r["ec50_2b"] for r in sel_results],
                    w3b, label="5-HT2B EC50 (nM)", color="#DD8452", alpha=0.85)
ax3b.set_yscale("log")
ax3b.set_xticks(x3b)
ax3b.set_xticklabels([c.replace("4C-D (Ariadne)","4C-D\n(Ariadne)")
                      for c in compound_names], fontsize=9)
ax3b.set_ylabel("EC50 (nM, log scale)")
ax3b.set_title("EC50 at 5-HT2A vs 5-HT2B\n(log scale)")
ax3b.legend(fontsize=8)

# Annotate fold selectivity
for xi, r in enumerate(sel_results):
    ax3b.text(xi, max(r["ec50_2a"], r["ec50_2b"]) * 1.3,
              f"{r['fold']:.0f}x", ha="center", fontsize=8, color="#333333")

fig3.suptitle("Ariadne 4-position analogs: 5-HT2A vs 5-HT2B selectivity (Gq BRET)\n"
              "Cardiac safety: 2A-selective = safer (5-HT2B agonism → valvulopathy)",
              fontsize=10)
plt.tight_layout()
plt.savefig("analysis/ariadne_2a2b_selectivity.png", dpi=150)
p("Saved analysis/ariadne_2a2b_selectivity.png")
plt.close()

log.close()
print("Done. Report: analysis/ariadne_report.txt")
print("Plots: ariadne_gq_vs_caflux.png | ariadne_alpha_series.png "
      "| ariadne_2a2b_selectivity.png")
