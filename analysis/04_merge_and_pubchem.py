"""
Script 04: Merge literature-curated data with ChEMBL core dataset.
Computes the proper log bias factor using both Emax and EC50.

Bias factor definition (operational model of agonism, Kenakin et al.):
  transduction ratio per pathway = Emax / EC50
  log_bias = log10(Emax_barr / EC50_barr) - log10(Emax_gq / EC50_gq)
           = log10( (Emax_barr / EC50_barr) / (Emax_gq / EC50_gq) )

Positive = barr-biased, Negative = Gq-biased, ~0 = balanced.
Threshold: |log_bias| > 0.5 (approx 3x preference) for classification.

Falls back to Emax-only ratio for compounds missing EC50 in one pathway.

Outputs:
  analysis/combined_dataset.csv
  analysis/04_merge_report.txt
  analysis/bias_combined.png
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── 1. Load literature data ───────────────────────────────────────────────────

lit = pd.read_csv("analysis/literature_curated.csv")
print(f"Loaded {len(lit)} literature compounds")
print(f"Columns: {list(lit.columns)}")

# ── 2. Compute bias factors for literature compounds ──────────────────────────

def compute_log_bias(row):
    """
    Prefer full transduction ratio (Emax/EC50) for both pathways.
    Fall back to Emax-only if EC50 missing for either pathway.
    Returns (log_bias, method_used).
    """
    gq_emax  = row.get("gq_emax")
    gq_ec50  = row.get("gq_ec50_nM")
    ba_emax  = row.get("barr_emax")
    ba_ec50  = row.get("barr_ec50_nM")

    # Full transduction ratio
    if (pd.notna(gq_emax) and pd.notna(gq_ec50) and
        pd.notna(ba_emax) and pd.notna(ba_ec50) and
        gq_emax > 0 and gq_ec50 > 0 and ba_ec50 > 0):
        tr_gq = gq_emax / gq_ec50
        tr_ba = ba_emax / ba_ec50
        return np.log10(tr_ba / tr_gq), "Emax/EC50"

    # Emax-only fallback
    if pd.notna(gq_emax) and pd.notna(ba_emax) and gq_emax > 0:
        return np.log10(ba_emax / gq_emax), "Emax_only"

    return np.nan, "no_data"

results = [compute_log_bias(row) for _, row in lit.iterrows()]
lit["log_bias"]    = [r[0] for r in results]
lit["bias_method"] = [r[1] for r in results]

def classify_bias(log_bias, threshold=0.5):
    if pd.isna(log_bias):   return "no data"
    if log_bias >  threshold: return "barr-biased"
    if log_bias < -threshold: return "Gq-biased"
    return "balanced"

lit["bias_class"] = lit["log_bias"].apply(classify_bias)
lit["source"] = "literature (Pottie & Poulie 2023)"

print("\nLiterature bias distribution (Emax/EC50 method):")
print(lit["bias_class"].value_counts().to_string())
print("\nLiterature log bias range:")
print(lit[["compound_name", "log_bias", "bias_class", "bias_method"]].sort_values("log_bias").to_string(index=False))

# ── 3. Load and recompute ChEMBL bias factors ─────────────────────────────────
# ChEMBL compounds only had Emax (EC50 was not saved to paired_compounds.csv).
# We reload from bias_factors.csv which has the Emax-based log_bias already.

chembl = pd.read_csv("analysis/bias_factors.csv")
chembl["source"]             = chembl["paper"]
chembl["compound_name"]      = chembl["pref_name"].fillna(chembl["molecule_chembl_id"])
chembl["reference_agonist"]  = "LSD/serotonin (per paper)"
chembl["compound_id_in_paper"] = chembl["molecule_chembl_id"]
chembl["bias_method"]        = "Emax_only"
# reclassify ChEMBL with same threshold
chembl["bias_class"] = chembl["log_bias"].apply(classify_bias)

print("\nChEMBL bias distribution (Emax-only, threshold=0.5):")
print(chembl["bias_class"].value_counts().to_string())

# ── 4. Align and merge ────────────────────────────────────────────────────────

shared_cols = ["compound_name", "compound_id_in_paper", "smiles",
               "gq_emax", "gq_ec50_nM", "barr_emax", "barr_ec50_nM",
               "log_bias", "bias_class", "bias_method",
               "reference_agonist", "source"]

chembl = chembl.rename(columns={"gq_ec50": "gq_ec50_nM", "barr_ec50": "barr_ec50_nM"})

def align(df, cols):
    for c in cols:
        if c not in df.columns:
            df[c] = None
    return df[cols].copy()

# Drop reference/control compounds from literature before merging
lit_clean = lit[~lit["compound_name"].isin({"LSD", "5-HT", "serotonin", "SEROTONIN"})]

lit_aligned    = align(lit_clean, shared_cols)
chembl_aligned = align(chembl,    shared_cols)

combined = pd.concat([chembl_aligned, lit_aligned], ignore_index=True)
combined = combined[combined["smiles"].notna()].reset_index(drop=True)

print(f"\nFinal combined dataset: {len(combined)} compounds")
print("\nBias class distribution:")
print(combined["bias_class"].value_counts().to_string())
print("\nBy source:")
print(combined.groupby("source")["bias_class"].value_counts().to_string())
print(f"\nLog bias  min={combined['log_bias'].min():.3f}  "
      f"max={combined['log_bias'].max():.3f}  "
      f"mean={combined['log_bias'].mean():.3f}  "
      f"std={combined['log_bias'].std():.3f}")

combined.to_csv("analysis/combined_dataset.csv", index=False)
print("\nSaved analysis/combined_dataset.csv")

# ── 5. Cross-check: compare our log_bias with paper's beta-factor ─────────────
# The paper reports beta for each compound. Let's see if our log_bias
# correlates with their beta (they should be monotonically related).

print("\n" + "="*70)
print("CROSS-CHECK: our log_bias vs paper beta-factor")
print("(beta column not in our CSV but we can compare qualitatively)")
print("Expected: higher beta -> more barr-biased -> higher log_bias")
print("="*70)
lit_check = lit_clean.sort_values("log_bias")
print(lit_check[["compound_name", "gq_emax", "gq_ec50_nM",
                  "barr_emax", "barr_ec50_nM",
                  "log_bias", "bias_class"]].to_string(index=False))

# ── 6. Plot ───────────────────────────────────────────────────────────────────

source_colors = {
    "25CN-NBOH (J Med Chem 2022)":                     "#4C72B0",
    "Heterocyclic Gq-biased (ACS Med Chem Lett 2023)":  "#DD8452",
    "literature (Pottie & Poulie 2023)":                "#55A868",
}

fig, axes = plt.subplots(1, 2, figsize=(13, 4))

ax = axes[0]
for source, grp in combined.groupby("source"):
    ax.hist(grp["log_bias"].dropna(), bins=10, alpha=0.7,
            label=source[:42], color=source_colors.get(source, "gray"))
ax.axvline(0,    color="black", linestyle="--", linewidth=1, label="balanced")
ax.axvline( 0.5, color="gray",  linestyle=":",  linewidth=1, label="±0.5 threshold")
ax.axvline(-0.5, color="gray",  linestyle=":",  linewidth=1)
ax.set_xlabel("log10 bias factor")
ax.set_ylabel("Count")
ax.set_title("Signalling bias — combined dataset\n(Emax/EC50 where available)")
ax.legend(fontsize=7)

ax = axes[1]
for source, grp in combined.groupby("source"):
    ax.scatter(grp["gq_emax"], grp["barr_emax"], alpha=0.8,
               label=source[:42], color=source_colors.get(source, "gray"))
lim = max(combined["gq_emax"].max(), combined["barr_emax"].max()) * 1.05
ax.plot([0, lim], [0, lim], "k--", linewidth=1, label="balanced (y=x)")
ax.set_xlabel("Gq Emax (%)")
ax.set_ylabel("β-arrestin Emax (%)")
ax.set_title("Gq vs β-arrestin efficacy — combined")
ax.legend(fontsize=7)

plt.tight_layout()
plt.savefig("analysis/bias_combined.png", dpi=150)
print("Saved analysis/bias_combined.png")

report = f"""
5-HT2A Biased Agonism — Combined Dataset Summary
=================================================
Bias factor: log10((Emax_barr/EC50_barr) / (Emax_gq/EC50_gq))
             Falls back to log10(Emax_barr/Emax_gq) when EC50 unavailable.
Classification threshold: |log_bias| > 0.5 (~3x pathway preference)

Sources:
  ChEMBL (2 papers, Emax-only):          {len(chembl_aligned)} compounds
  Literature (Pottie & Poulie 2023):     {len(lit_aligned)} compounds

Final combined dataset:                  {len(combined)} compounds

Bias class distribution:
{combined['bias_class'].value_counts().to_string()}

By source:
{combined.groupby('source')['bias_class'].value_counts().to_string()}

Log bias statistics:
  Min:  {combined['log_bias'].min():.3f}
  Max:  {combined['log_bias'].max():.3f}
  Mean: {combined['log_bias'].mean():.3f}
  Std:  {combined['log_bias'].std():.3f}

Bias method per source:
  Pottie & Poulie 2023: Emax/EC50 transduction ratio
  ChEMBL compounds:     Emax ratio only (EC50 not available in paired_compounds.csv)

Note: The ChEMBL compounds should ideally be re-fetched with EC50 data
to use the full transduction ratio. This is a next-step improvement.

Next steps:
  - Compute Morgan fingerprints for all 46 compounds
  - Build QSAR model predicting log_bias from fingerprints
  - Evaluate with leave-one-out or 5-fold cross-validation
"""
print(report)
with open("analysis/04_merge_report.txt", "w") as f:
    f.write(report)
print("Done.")
