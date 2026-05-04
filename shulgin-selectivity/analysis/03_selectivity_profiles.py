"""
Script 03: Compute serotonin receptor selectivity profiles for Shulgin compounds.

For each compound with data at >= 3 receptors:
  1. Compute geometric mean pKi per receptor (across replicates/species)
  2. Apply four selectivity definitions from the kinase selectivity framework:
       - S-score  : fraction of receptors below a potency threshold
       - Gini     : concentration of binding in the top receptor
       - Entropy  : Shannon entropy of the pKi distribution
       - Ratio    : pKi_top1 / pKi_top2
  3. Compute 5-HT2A selectivity over 5-HT2B specifically (safety-relevant)

Outputs:
  analysis/selectivity_profiles.csv    - one row per compound
  analysis/pki_matrix.csv             - pKi values per compound x receptor
  analysis/03_selectivity_report.txt  - plain-text summary
  analysis/selectivity_heatmap.png    - heatmap of pKi profiles
  analysis/selectivity_scores.png     - scatter plots of selectivity metrics
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

log = open("analysis/03_selectivity_report.txt", "w")
def p(*args, **kwargs):
    print(*args, **kwargs, file=log)

# ── 1. Load clean data ────────────────────────────────────────────────────────

df = pd.read_csv("analysis/shulgin_clean.csv", low_memory=False)
df = df[df["ki_nM"].notna() & ~df["censored"]].copy()
p(f"Loaded {len(df)} valid (non-censored) Ki records")
p(f"Compounds: {df['canonical_name'].nunique()}")

RECEPTORS = ["5-HT1A", "5-HT1B", "5-HT1D", "5-HT1E", "5-HT1F",
             "5-HT2A", "5-HT2B", "5-HT2C", "5-HT3", "5-HT4",
             "5-HT6", "5-HT7"]

# ── 2. Geometric mean pKi per compound x receptor ────────────────────────────
# Geometric mean in pKi space = arithmetic mean of pKi values
# (equivalent to geometric mean of Ki values)

pki_rows = []
for compound, grp in df.groupby("canonical_name"):
    row = {"compound": compound}
    for rec in RECEPTORS:
        rec_data = grp[grp["name"] == rec]["pKi"]
        if len(rec_data) > 0:
            row[rec] = rec_data.mean()  # mean of pKi = geometric mean of Ki
        else:
            row[rec] = np.nan
    pki_rows.append(row)

pki_df = pd.DataFrame(pki_rows).set_index("compound")
pki_df.to_csv("analysis/pki_matrix.csv")
p(f"\npKi matrix: {pki_df.shape[0]} compounds x {pki_df.shape[1]} receptors")

# Filter to compounds with >= 3 receptors
n_receptors = pki_df.notna().sum(axis=1)
pki_df = pki_df[n_receptors >= 3]
p(f"After filtering (>= 3 receptors): {len(pki_df)} compounds")

# ── 3. Selectivity definitions ────────────────────────────────────────────────

THRESHOLD = 6.0  # pKi threshold for S-score (~1 uM)

def s_score(row):
    """Fraction of tested receptors BELOW threshold (higher = more selective)."""
    vals = row.dropna()
    if len(vals) == 0: return np.nan
    return (vals < THRESHOLD).sum() / len(vals)

def gini(row):
    """Gini coefficient of pKi distribution (higher = more concentrated/selective)."""
    vals = row.dropna().values
    if len(vals) < 2: return np.nan
    vals = np.sort(vals)
    n = len(vals)
    idx = np.arange(1, n + 1)
    return (2 * np.sum(idx * vals) / (n * np.sum(vals))) - (n + 1) / n

def entropy(row):
    """Shannon entropy of normalised pKi distribution (lower = more selective)."""
    vals = row.dropna().values
    if len(vals) < 2: return np.nan
    vals = vals - vals.min() + 1e-9  # shift to positive
    probs = vals / vals.sum()
    return -np.sum(probs * np.log2(probs + 1e-12))

def ratio(row):
    """Ratio of top pKi to second-highest pKi (higher = more selective)."""
    vals = row.dropna().sort_values(ascending=False)
    if len(vals) < 2: return np.nan
    if vals.iloc[1] == 0: return np.nan
    return vals.iloc[0] / vals.iloc[1]

scores = pd.DataFrame({
    "n_receptors":  pki_df.notna().sum(axis=1),
    "s_score":      pki_df.apply(s_score,  axis=1),
    "gini":         pki_df.apply(gini,     axis=1),
    "entropy":      pki_df.apply(entropy,  axis=1),
    "ratio":        pki_df.apply(ratio,    axis=1),
    "top_receptor": pki_df.apply(lambda r: r.dropna().idxmax() if r.notna().any() else np.nan, axis=1),
    "pKi_top":      pki_df.apply(lambda r: r.dropna().max()    if r.notna().any() else np.nan, axis=1),
    "pKi_2A":       pki_df["5-HT2A"],
    "pKi_2B":       pki_df["5-HT2B"],
})

# 5-HT2A vs 5-HT2B selectivity (cardiac safety)
scores["selectivity_2A_over_2B"] = scores["pKi_2A"] - scores["pKi_2B"]
# Positive = 2A selective (safer), negative = 2B selective (cardiac risk)

# ── 4. Report ─────────────────────────────────────────────────────────────────

p(f"\n{'='*75}")
p(f"SELECTIVITY PROFILES — Shulgin compounds at serotonin receptors")
p(f"{'='*75}")
p(f"\nDefinitions:")
p(f"  S-score:  fraction of receptors below pKi={THRESHOLD} (higher = more selective)")
p(f"  Gini:     concentration of binding (higher = more selective)")
p(f"  Entropy:  Shannon entropy of pKi distribution (lower = more selective)")
p(f"  Ratio:    pKi_top1 / pKi_top2 (higher = more selective)")
p(f"  2A-2B:    pKi(5-HT2A) - pKi(5-HT2B) (positive = 2A-selective = safer)")

p(f"\n{'Compound':<22} {'nRec':>4} {'S-sc':>6} {'Gini':>6} {'Entr':>6} "
  f"{'Ratio':>6} {'2A-2B':>6}  Top receptor")
p("-"*85)

for compound, row in scores.sort_values("selectivity_2A_over_2B",
                                         ascending=False).iterrows():
    p(f"  {compound:<20} {row['n_receptors']:>4.0f} "
      f"{row['s_score']:>6.3f} {row['gini']:>6.3f} "
      f"{row['entropy']:>6.3f} {row['ratio']:>6.3f} "
      f"{row['selectivity_2A_over_2B']:>+6.2f}  {row['top_receptor']}")

# Most/least 2A-selective
p(f"\n{'='*50}")
p(f"5-HT2A vs 5-HT2B selectivity (cardiac safety)")
p(f"{'='*50}")
p(f"Most 2A-selective (safest):")
top5 = scores.dropna(subset=["selectivity_2A_over_2B"]).nlargest(5, "selectivity_2A_over_2B")
for c, r in top5.iterrows():
    p(f"  {c:<22} 2A-2B = {r['selectivity_2A_over_2B']:+.2f}  "
      f"pKi_2A={r['pKi_2A']:.2f}  pKi_2B={r['pKi_2B']:.2f}")
p(f"\nLeast 2A-selective (most cardiac risk):")
bot5 = scores.dropna(subset=["selectivity_2A_over_2B"]).nsmallest(5, "selectivity_2A_over_2B")
for c, r in bot5.iterrows():
    p(f"  {c:<22} 2A-2B = {r['selectivity_2A_over_2B']:+.2f}  "
      f"pKi_2A={r['pKi_2A']:.2f}  pKi_2B={r['pKi_2B']:.2f}")

# Correlation between definitions
p(f"\n{'='*50}")
p(f"Correlation between selectivity definitions")
p(f"{'='*50}")
corr_cols = ["s_score", "gini", "entropy", "ratio"]
corr = scores[corr_cols].corr()
p(corr.round(3).to_string())

# ── 5. Save scores ────────────────────────────────────────────────────────────

full = pki_df.copy()
full.index.name = "compound"
full = full.merge(scores, left_index=True, right_index=True)
full.to_csv("analysis/selectivity_profiles.csv")
p(f"\nSaved analysis/selectivity_profiles.csv")

# ── 6. Heatmap ────────────────────────────────────────────────────────────────

# Only compounds with >= 5 receptors for the heatmap
heatmap_df = pki_df[pki_df.notna().sum(axis=1) >= 5]
# Sort by 5-HT2A pKi
heatmap_df = heatmap_df.reindex(
    heatmap_df["5-HT2A"].sort_values(ascending=False).index)

fig, ax = plt.subplots(figsize=(13, max(6, len(heatmap_df) * 0.4)))
im = ax.imshow(heatmap_df.values, aspect="auto", cmap="YlOrRd",
               vmin=4, vmax=10)
ax.set_xticks(range(len(heatmap_df.columns)))
ax.set_xticklabels(heatmap_df.columns, rotation=45, ha="right", fontsize=9)
ax.set_yticks(range(len(heatmap_df)))
ax.set_yticklabels(heatmap_df.index, fontsize=8)

# Annotate cells
for i in range(len(heatmap_df)):
    for j in range(len(heatmap_df.columns)):
        val = heatmap_df.values[i, j]
        if not np.isnan(val):
            ax.text(j, i, f"{val:.1f}", ha="center", va="center",
                    fontsize=6, color="black" if val < 8 else "white")

plt.colorbar(im, ax=ax, label="pKi")
ax.set_title("Serotonin receptor binding profiles — Shulgin compounds\n"
             "(pKi = -log10[Ki/nM] + 9; higher = more potent; grey = no data)",
             fontsize=10)
plt.tight_layout()
plt.savefig("analysis/selectivity_heatmap.png", dpi=150)
p("Saved analysis/selectivity_heatmap.png")

# ── 7. Selectivity score scatter plots ────────────────────────────────────────

fig2, axes = plt.subplots(1, 3, figsize=(14, 5))

pairs = [("s_score", "gini"), ("s_score", "entropy"), ("gini", "entropy")]
for ax, (x, y) in zip(axes, pairs):
    sc = scores.dropna(subset=[x, y])
    ax.scatter(sc[x], sc[y], alpha=0.7, s=50, color="#4C72B0")
    for name, row in sc.iterrows():
        ax.annotate(name[:12], (row[x], row[y]),
                    fontsize=5, alpha=0.7,
                    xytext=(3, 3), textcoords="offset points")
    r = sc[[x, y]].corr().iloc[0, 1]
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.set_title(f"{x} vs {y}\nr = {r:.3f}")

plt.suptitle("Selectivity definition correlations — Shulgin compounds", fontsize=10)
plt.tight_layout()
plt.savefig("analysis/selectivity_scores.png", dpi=150)
p("Saved analysis/selectivity_scores.png")

log.close()
print("Done. Report: analysis/03_selectivity_report.txt")
print(f"Compounds analysed: {len(scores)}")
print(f"Plots: selectivity_heatmap.png, selectivity_scores.png")
