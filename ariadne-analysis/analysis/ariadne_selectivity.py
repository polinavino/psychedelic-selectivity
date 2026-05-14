"""
Selectivity profile analysis of Ariadne and the alpha-alkyl DOx series.

Data source: Cunningham et al. ACS Chem Neurosci 2023
             doi:10.1021/acschemneuro.2c00597
             Loaded from: data/cunningham2023_data.csv

Parameters:
  SSCORE_THRESHOLD : pEC50 below which a receptor counts as inactive for S-score
  PRIMARY_ASSAY    : assay used for selectivity profiles ("gq_bret" or "ca_flux")

Outputs:
  analysis/ariadne_selectivity_report.txt
  analysis/ariadne_gq_vs_caflux.png      - Gq vs Ca2+ flux: potency + efficacy
  analysis/ariadne_alpha_series.png       - alpha-alkyl series across 4-position pairs
  analysis/ariadne_profiles.png           - overview pEC50 + alpha-alkyl trend
  analysis/ariadne_selectivity_scores.png - four selectivity metrics
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import os

os.makedirs("analysis", exist_ok=True)

# ── Parameters ────────────────────────────────────────────────────────────────

SSCORE_THRESHOLD = 6.5
PRIMARY_ASSAY    = "gq_bret"

log = open("analysis/ariadne_selectivity_report.txt", "w")
def p(*args, **kwargs):
    print(*args, **kwargs, file=log)

# ── 1. Load data ──────────────────────────────────────────────────────────────

df = pd.read_csv("data/cunningham2023_data.csv", comment="#")
df["ec50_nM"]   = pd.to_numeric(df["ec50_nM"],  errors="coerce")
df["emax_pct"]  = pd.to_numeric(df["emax_pct"], errors="coerce")
df["pec50"]     = 9 - np.log10(df["ec50_nM"])
df["log_trans"] = np.log10(df["emax_pct"] / df["ec50_nM"])

p("Ariadne Selectivity Analysis — Cunningham et al. ACS Chem Neurosci 2023")
p("=" * 70)
p(f"Data file:       data/cunningham2023_data.csv")
p(f"Primary assay:   {PRIMARY_ASSAY}")
p(f"S-score threshold: pEC50 < {SSCORE_THRESHOLD} "
  f"(EC50 > {10**(9-SSCORE_THRESHOLD):.1f} nM)")
p(f"Records loaded:  {len(df)}")

# ── 2. Selectivity metric functions ──────────────────────────────────────────

def s_score(vals):
    vals = [v for v in vals if pd.notna(v)]
    if not vals: return np.nan
    return sum(v < SSCORE_THRESHOLD for v in vals) / len(vals)

def gini(vals):
    vals = sorted([v for v in vals if pd.notna(v)])
    if len(vals) < 2: return np.nan
    n = len(vals); idx = np.arange(1, n+1)
    return (2 * np.sum(idx * np.array(vals)) / (n * np.sum(vals))) - (n+1)/n

def entropy(vals):
    vals = np.array([v for v in vals if pd.notna(v)])
    if len(vals) < 2: return np.nan
    vals = vals - vals.min() + 1e-9
    probs = vals / vals.sum()
    return -np.sum(probs * np.log2(probs + 1e-12))

def ratio(vals):
    vals = sorted([v for v in vals if pd.notna(v)], reverse=True)
    if len(vals) < 2 or vals[1] == 0: return np.nan
    return vals[0] / vals[1]

# ── 3. Selectivity profiles ───────────────────────────────────────────────────

profile_compounds = ["5-HT", "(rac)-Ariadne", "(R)-Ariadne", "(S)-Ariadne"]
receptors_2ht     = ["5-HT2A", "5-HT2B", "5-HT2C"]

p(f"\n\n{'='*70}")
p(f"Selectivity profiles — {PRIMARY_ASSAY}, 5-HT2A/2B/2C receptors")
p("=" * 70)
p(f"\n{'Compound':<20} {'S-score':>8} {'Gini':>8} {'Entropy':>8} "
  f"{'Ratio':>8}  Top receptor")
p("-" * 65)

profile_scores = []
for compound in profile_compounds:
    sub  = df[(df["compound"]==compound) & (df["assay"]==PRIMARY_ASSAY) &
              (df["receptor"].isin(receptors_2ht))]
    vals = sub["pec50"].tolist()
    ss   = s_score(vals)
    gi   = gini(vals)
    en   = entropy(vals)
    ra   = ratio(vals)
    top  = sub.loc[sub["pec50"].idxmax(), "receptor"] if len(sub) > 0 else "n/a"
    p(f"  {compound:<18} {ss:>8.3f} {gi:>8.3f} {en:>8.3f} {ra:>8.3f}  {top}")
    profile_scores.append({"compound": compound, "s_score": ss, "gini": gi,
                            "entropy": en, "ratio": ra, "top_receptor": top})

sc_df = pd.DataFrame(profile_scores)

p(f"\nCorrelations between selectivity definitions (n=4 — direction only):")
p(sc_df[["s_score","gini","entropy","ratio"]].corr().round(3).to_string())

# ── 4. Alpha-alkyl series ─────────────────────────────────────────────────────

p(f"\n\n{'='*70}")
p("Alpha-alkyl series at 5-HT2A (Gq and beta-arrestin2)")
p("=" * 70)
p(f"\n{'Compound':<25} {'alpha':>5} {'Gq EC50':>9} {'Gq Emax':>9} "
  f"{'barr EC50':>10} {'barr Emax':>10} {'logBias':>9} {'Halluc?':>8}")
p("-" * 92)

alpha_series = {
    "DOM":           ("Me",  True),
    "(rac)-Ariadne": ("Et",  False),
    "5C-D":          ("nPr", False),
    "4C-Pr":         ("Et",  False),
    "DOPr":          ("Me",  True),
    "4C-I":          ("Et",  False),
    "DOI":           ("Me",  True),
}

for compound, (alpha, halluc) in alpha_series.items():
    gq   = df[(df["compound"]==compound) & (df["assay"]=="gq_bret") &
              (df["receptor"]=="5-HT2A")]
    barr = df[(df["compound"]==compound) & (df["assay"]=="barr2_bret") &
              (df["receptor"]=="5-HT2A")]
    gq_ec50 = gq["ec50_nM"].values[0]   if len(gq)   > 0 else np.nan
    gq_emax = gq["emax_pct"].values[0]  if len(gq)   > 0 else np.nan
    ba_ec50 = barr["ec50_nM"].values[0] if len(barr) > 0 else np.nan
    ba_emax = barr["emax_pct"].values[0]if len(barr) > 0 else np.nan
    if pd.notna(gq_ec50) and pd.notna(ba_ec50):
        lb   = np.log10((ba_emax/ba_ec50) / (gq_emax/gq_ec50))
        lb_s = f"{lb:+9.3f}"
    else:
        lb_s = "      n/a"
    p(f"  {compound:<25} {alpha:>5} {gq_ec50:>9.1f} {gq_emax:>9.1f} "
      f"{ba_ec50 if pd.notna(ba_ec50) else float('nan'):>10.1f} "
      f"{ba_emax if pd.notna(ba_emax) else float('nan'):>10.1f} "
      f"{lb_s} {'YES' if halluc else 'no':>8}")

# ── 5. 4-position analog selectivity ─────────────────────────────────────────

p(f"\n\n{'='*70}")
p("4-position analogs: 5-HT2A vs 5-HT2B selectivity (Gq BRET)")
p("=" * 70)
p(f"\n{'Compound':<20} {'2A EC50':>9} {'2A Emax':>9} {'2B EC50':>9} "
  f"{'2B Emax':>9} {'pEC50 2A-2B':>12}")
p("-" * 72)

for compound in ["4C-D (Ariadne)","4C-TFM","4C-Pr","4C-cycPr","4C-MOM"]:
    s2a = df[(df["compound"]==compound)&(df["assay"]=="gq_bret")&
             (df["receptor"]=="5-HT2A")]
    s2b = df[(df["compound"]==compound)&(df["assay"]=="gq_bret")&
             (df["receptor"]=="5-HT2B")]
    if len(s2a)==0 or len(s2b)==0: continue
    e2a = s2a["ec50_nM"].values[0]; m2a = s2a["emax_pct"].values[0]
    e2b = s2b["ec50_nM"].values[0]; m2b = s2b["emax_pct"].values[0]
    sel = (9-np.log10(e2a)) - (9-np.log10(e2b))
    p(f"  {compound:<20} {e2a:>9.1f} {m2a:>9.1f} {e2b:>9.1f} "
      f"{m2b:>9.1f} {sel:>+12.2f}")

p(f"\n  4C-TFM: 8-fold 2A/2B selectivity "
  f"(pEC50 diff = {np.log10(234.0/28.6):.2f}) — most selective analog.")

# ── Colour map ────────────────────────────────────────────────────────────────

colors = {
    "5-HT":          "#888888",
    "(rac)-Ariadne": "#4C72B0",
    "(R)-Ariadne":   "#2196F3",
    "(S)-Ariadne":   "#90CAF9",
    "DOM":           "#DD8452",
}

# ── Plot 1: Gq BRET vs Ca2+ flux — 2x3 grid ──────────────────────────────────

fig1, axes1 = plt.subplots(2, 3, figsize=(14, 8), sharey="row")
assay_labels = {"gq_bret": "Gq dissociation BRET", "ca_flux": "Ca²⁺ flux"}

for row, assay in enumerate(["gq_bret", "ca_flux"]):
    for col, receptor in enumerate(receptors_2ht):
        ax  = axes1[row][col]
        sub = df[(df["assay"]==assay) & (df["receptor"]==receptor) &
                 (df["compound"].isin(profile_compounds))]
        sub = sub.set_index("compound").reindex(profile_compounds)
        x   = np.arange(len(profile_compounds))
        w   = 0.5
        bar_colors = [colors.get(c, "gray") for c in profile_compounds]

        ax.bar(x, sub["pec50"].values, w, color=bar_colors, alpha=0.8)

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
        ax.axhline(SSCORE_THRESHOLD, color="red", lw=0.8, linestyle="--", alpha=0.5)

        for xi, compound in enumerate(profile_compounds):
            emax = sub.loc[compound, "emax_pct"] if compound in sub.index else np.nan
            if pd.notna(emax) and compound != "5-HT":
                ax2.text(xi, emax+4, f"{emax:.0f}%",
                         ha="center", fontsize=6, color="black")

legend_els = [Patch(color=colors[c], alpha=0.8, label=c)
              for c in profile_compounds]
fig1.legend(handles=legend_els, loc="lower center", ncol=4,
            fontsize=8, bbox_to_anchor=(0.5, -0.01))
fig1.suptitle("Ariadne: potency (pEC50, bars) and efficacy (Emax %, line)\n"
              "Gq BRET (top row) vs Ca²⁺ flux (bottom row)",
              fontsize=10)
plt.tight_layout(rect=[0, 0.05, 1, 1])
plt.savefig("analysis/ariadne_gq_vs_caflux.png", dpi=150, bbox_inches="tight")
p("Saved analysis/ariadne_gq_vs_caflux.png")
plt.close()

# ── Plot 2: alpha-alkyl bias series — 3 panels ───────────────────────────────

pairs = [
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

for ax, (title, members) in zip(axes2, pairs):
    labels2, pec50s, emaxs, log_biases, bar_cols = [], [], [], [], []
    for compound, alpha, halluc, color in members:
        gq_row = df[(df["compound"]==compound) & (df["assay"]=="gq_bret") &
                    (df["receptor"]=="5-HT2A")]
        ba_row = df[(df["compound"]==compound) & (df["assay"]=="barr2_bret") &
                    (df["receptor"]=="5-HT2A")]
        if len(gq_row) == 0: continue
        ec50 = gq_row["ec50_nM"].values[0]
        emax = gq_row["emax_pct"].values[0]
        pec50s.append(9 - np.log10(ec50))
        emaxs.append(emax)
        bar_cols.append(color)
        if len(ba_row) > 0:
            bec50 = ba_row["ec50_nM"].values[0]
            bemax = ba_row["emax_pct"].values[0]
            log_biases.append(np.log10((bemax/bec50) / (emax/ec50)))
        else:
            log_biases.append(np.nan)
        hmark = "★" if halluc else "○"
        labels2.append(f"{hmark} {compound}\n(α-{alpha})")

    x    = np.arange(len(labels2))
    bars = ax.bar(x, pec50s, color=bar_cols, alpha=0.85, width=0.5)
    ax_t = ax.twinx()
    ax_t.plot(x, emaxs, "ko--", markersize=7, linewidth=1.5)
    ax_t.set_ylim(0, 120)
    ax_t.set_ylabel("Emax (% 5-HT)", fontsize=8)

    for xi, (lb, pe) in enumerate(zip(log_biases, pec50s)):
        if pd.notna(lb):
            ax.text(xi, pe+0.15, f"bias={lb:+.2f}",
                    ha="center", va="bottom", fontsize=7.5, color="#333333")

    ax.set_xticks(x)
    ax.set_xticklabels(labels2, fontsize=8)
    ax.set_ylabel("pEC50 at 5-HT2A (Gq)", fontsize=8)
    ax.set_ylim(0, 10)
    ax.set_title(f"{title}\nbars=pEC50  line=Emax  label=log_bias", fontsize=8)

fig2.suptitle("Alpha-substituent effect across 4-position pairs\n"
              "★ = hallucinogenic (alpha-methyl)   ○ = non-hallucinogenic (alpha-ethyl)\n"
              "log_bias: negative = Gq-biased; note consistent small shift, not reversal",
              fontsize=9)
plt.tight_layout()
plt.savefig("analysis/ariadne_alpha_series.png", dpi=150)
p("Saved analysis/ariadne_alpha_series.png")
plt.close()

# ── Plot 3: ariadne_profiles.png — overview ───────────────────────────────────

fig3, axes3 = plt.subplots(1, 2, figsize=(13, 5))
ax = axes3[0]
x  = np.arange(len(receptors_2ht))
w  = 0.18
for i, compound in enumerate(profile_compounds):
    sub  = df[(df["compound"]==compound) & (df["assay"]==PRIMARY_ASSAY) &
              (df["receptor"].isin(receptors_2ht))].set_index("receptor")
    vals = [sub.loc[r,"pec50"] if r in sub.index else np.nan for r in receptors_2ht]
    ax.bar(x+(i-1.5)*w, vals, w,
           label=compound, color=colors.get(compound,"gray"), alpha=0.85)
ax.axhline(SSCORE_THRESHOLD, color="red", linestyle="--", lw=1,
           label=f"S-score threshold (pEC50={SSCORE_THRESHOLD})")
ax.set_xticks(x); ax.set_xticklabels(receptors_2ht)
ax.set_ylabel("pEC50 (Gq dissociation BRET)")
ax.set_title("5-HT2A/2B/2C potency: Ariadne vs 5-HT")
ax.legend(fontsize=7); ax.set_ylim(0, 12)

ax2 = axes3[1]
series3 = [("DOM","Me",True,"#DD8452"),
           ("(rac)-Ariadne","Et",False,"#4C72B0"),
           ("5C-D","nPr",False,"#55A868")]
pe3,em3,lb3,c3 = [],[],[],[]
for compound,alpha,halluc,color in series3:
    sub = df[(df["compound"]==compound)&(df["assay"]=="gq_bret")&
             (df["receptor"]=="5-HT2A")]
    if len(sub)==0: continue
    pe3.append(sub["pec50"].values[0])
    em3.append(sub["emax_pct"].values[0])
    lb3.append(f"{compound}\n(α-{alpha})\n{'Halluc.' if halluc else 'Non-halluc.'}")
    c3.append(color)
x3   = np.arange(len(lb3))
bars3= ax2.bar(x3, pe3, color=c3, alpha=0.85, width=0.5)
ax2t = ax2.twinx()
ax2t.plot(x3, em3, "ko--", markersize=8)
ax2t.set_ylabel("Emax (% 5-HT)"); ax2t.set_ylim(0,120)
for bar,v in zip(bars3,pe3):
    ax2.text(bar.get_x()+bar.get_width()/2.,v+0.05,
             f"{v:.2f}",ha="center",va="bottom",fontsize=9)
ax2.set_xticks(x3); ax2.set_xticklabels(lb3, fontsize=8)
ax2.set_ylabel("pEC50 at 5-HT2A (Gq)")
ax2.set_title("Alpha-chain length: potency & efficacy trend")
ax2.set_ylim(0,12)
plt.suptitle("Ariadne series overview — Cunningham et al. 2023", fontsize=10)
plt.tight_layout()
plt.savefig("analysis/ariadne_profiles.png", dpi=150)
p("Saved analysis/ariadne_profiles.png")
plt.close()

# ── Plot 4: selectivity scores ────────────────────────────────────────────────

fig4, axes4 = plt.subplots(1, 4, figsize=(14, 4))
for ax, sname, slabel in zip(axes4,
    ["s_score","gini","entropy","ratio"],
    ["S-score","Gini","Entropy (bits)","Ratio top1/top2"]):
    bar_colors = [colors.get(c,"gray") for c in sc_df["compound"]]
    ax.bar(range(len(sc_df)), sc_df[sname].values, color=bar_colors, alpha=0.85)
    ax.set_xticks(range(len(sc_df)))
    ax.set_xticklabels([c.replace("(","").replace(")","")
                        for c in sc_df["compound"]],
                       rotation=30, ha="right", fontsize=7)
    ax.set_title(slabel); ax.set_ylabel(slabel)
plt.suptitle("Selectivity definitions: Ariadne vs 5-HT\n"
             "(3 receptors: 5-HT2A/2B/2C — Gq BRET)", fontsize=9)
plt.tight_layout()
plt.savefig("analysis/ariadne_selectivity_scores.png", dpi=150)
p("Saved analysis/ariadne_selectivity_scores.png")
plt.close()

log.close()
print("Done. Outputs in analysis/")
print("  ariadne_selectivity_report.txt")
print("  ariadne_gq_vs_caflux.png")
print("  ariadne_alpha_series.png")
print("  ariadne_profiles.png")
print("  ariadne_selectivity_scores.png")
