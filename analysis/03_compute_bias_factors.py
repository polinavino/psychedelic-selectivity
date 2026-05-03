"""
Script 03: Compute bias factors for the 29 paired ChEMBL compounds.

For each compound we compute:
  - emax_ratio = Emax_barr / Emax_Gq
  - log_bias   = log10(Emax_barr / Emax_Gq)
    > 0  : beta-arrestin biased
    < 0  : Gq biased
    ~ 0  : balanced

EC50 is not used in the bias factor because the two source papers use
different reference agonists, making transduction ratios incomparable
across papers without re-normalisation.

Outputs:
  analysis/bias_factors.csv
  analysis/03_bias_report.txt
  analysis/bias_distribution.png
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Load paired records ───────────────────────────────────────────────────────

df_all = pd.read_csv("analysis/paired_compounds.csv")
print(f"Loaded {len(df_all)} records")
print(f"Standard types: {df_all['standard_type'].unique()}")
print(f"Pathways: {df_all['pathway'].unique()}")

names = pd.read_csv("analysis/paired_compound_names.csv")[
    ["molecule_chembl_id", "pref_name", "smiles", "mw", "alogp"]
]

# ── Extract Emax per compound per pathway ─────────────────────────────────────

def mean_by_compound(df, pathway, stype):
    subset = df[(df["pathway"] == pathway) & (df["standard_type"] == stype)]
    return subset.groupby("molecule_chembl_id")["standard_value"].mean()

gq_emax   = mean_by_compound(df_all, "gq",            "EMAX")
barr_emax = mean_by_compound(df_all, "beta-arrestin", "EMAX")

print(f"\nCompounds with Gq Emax:   {len(gq_emax)}")
print(f"Compounds with barr Emax: {len(barr_emax)}")

doc_map = (df_all[["molecule_chembl_id", "document_chembl_id"]]
           .drop_duplicates("molecule_chembl_id")
           .set_index("molecule_chembl_id")["document_chembl_id"])

result = pd.DataFrame({
    "gq_emax":   gq_emax,
    "barr_emax": barr_emax,
    "document_chembl_id": doc_map,
})
result.index.name = "molecule_chembl_id"
result = result.reset_index()

# ── Compute bias factors ──────────────────────────────────────────────────────

mask = result["gq_emax"].notna() & result["barr_emax"].notna() & (result["gq_emax"] > 0)
result.loc[mask, "emax_ratio"] = result.loc[mask, "barr_emax"] / result.loc[mask, "gq_emax"]
result.loc[mask, "log_bias"]   = np.log10(result.loc[mask, "emax_ratio"])

def classify_bias(log_bias):
    if pd.isna(log_bias): return "no data"
    if log_bias >  0.3:   return "barr-biased"
    if log_bias < -0.3:   return "Gq-biased"
    return "balanced"

result["bias_class"] = result["log_bias"].apply(classify_bias)

# ── Merge names and paper labels ──────────────────────────────────────────────

result = result.merge(names, on="molecule_chembl_id", how="left")

paper_labels = {
    "CHEMBL5113537": "25CN-NBOH (J Med Chem 2022)",
    "CHEMBL5370684": "Heterocyclic Gq-biased (ACS Med Chem Lett 2023)",
}
result["paper"] = result["document_chembl_id"].map(paper_labels)

# ── Print table ───────────────────────────────────────────────────────────────

def display_name(row):
    n = row["pref_name"]
    return str(n)[:21] if pd.notna(n) else row["molecule_chembl_id"][:21]

print("\n" + "="*85)
print(f"{'Compound':<22} {'Paper':<32} {'Gq Emax':>8} {'bArr Emax':>10} {'logBias':>8}  Class")
print("="*85)
for _, row in result.sort_values("log_bias").iterrows():
    print(f"{display_name(row):<22} {str(row['paper']):<32} "
          f"{row['gq_emax']:>8.1f} {row['barr_emax']:>10.1f} "
          f"{row['log_bias']:>8.3f}  {row['bias_class']}")

print("\nBias class distribution:")
print(result["bias_class"].value_counts().to_string())
print("\nBy paper:")
print(result.groupby("paper")["bias_class"].value_counts().to_string())
print(f"\nLog bias  min={result['log_bias'].min():.3f}  "
      f"max={result['log_bias'].max():.3f}  "
      f"mean={result['log_bias'].mean():.3f}  "
      f"std={result['log_bias'].std():.3f}")

# ── Plot ──────────────────────────────────────────────────────────────────────

colors = {
    "25CN-NBOH (J Med Chem 2022)":                    "#4C72B0",
    "Heterocyclic Gq-biased (ACS Med Chem Lett 2023)": "#DD8452",
}

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

ax = axes[0]
for paper, grp in result.groupby("paper"):
    ax.hist(grp["log_bias"].dropna(), bins=8, alpha=0.7,
            label=paper, color=colors.get(paper, "gray"))
ax.axvline(0,    color="black", linestyle="--", linewidth=1, label="balanced")
ax.axvline( 0.3, color="gray",  linestyle=":",  linewidth=1)
ax.axvline(-0.3, color="gray",  linestyle=":",  linewidth=1)
ax.set_xlabel("log10(Emax_βarr / Emax_Gq)")
ax.set_ylabel("Count")
ax.set_title("Signalling bias distribution")
ax.legend(fontsize=8)

ax = axes[1]
for paper, grp in result.groupby("paper"):
    ax.scatter(grp["gq_emax"], grp["barr_emax"], alpha=0.8,
               label=paper, color=colors.get(paper, "gray"))
lim = max(result["gq_emax"].max(), result["barr_emax"].max()) * 1.05
ax.plot([0, lim], [0, lim], "k--", linewidth=1, label="balanced (y=x)")
ax.set_xlabel("Gq Emax (%)")
ax.set_ylabel("β-arrestin Emax (%)")
ax.set_title("Gq vs β-arrestin efficacy")
ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig("analysis/bias_distribution.png", dpi=150)
print("\nSaved analysis/bias_distribution.png")

# ── Save ──────────────────────────────────────────────────────────────────────

result.to_csv("analysis/bias_factors.csv", index=False)
print("Saved analysis/bias_factors.csv")

report = f"""
5-HT2A Bias Factor Summary - ChEMBL Core Dataset
=================================================

Compounds with computable bias factor: {result['log_bias'].notna().sum()} / {len(result)}

Bias direction (threshold |log10 ratio| > 0.3, i.e. >2x preference):
{result['bias_class'].value_counts().to_string()}

Log bias statistics:
  Min:  {result['log_bias'].min():.3f}
  Max:  {result['log_bias'].max():.3f}
  Mean: {result['log_bias'].mean():.3f}
  Std:  {result['log_bias'].std():.3f}

By paper:
{result.groupby('paper')['bias_class'].value_counts().to_string()}

Note on EC50:
  The two source papers use different reference agonists (LSD vs serotonin),
  making EC50-based transduction ratios incomparable across papers without
  re-normalisation. Emax ratio is used as a more robust alternative.

Next steps:
  - Supplement with literature data to reach ~80-100 compounds
  - Compute Morgan fingerprints for all compounds
  - Build QSAR model predicting log_bias from molecular structure
"""
print(report)
with open("analysis/03_bias_report.txt", "w") as f:
    f.write(report)
print("Done.")
