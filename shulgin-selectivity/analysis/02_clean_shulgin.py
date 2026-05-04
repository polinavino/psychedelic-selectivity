"""
Script 02: Clean PDSP hits - keep only genuine Shulgin compounds.

The substring search in script 01 produced false positives (RISPERIDONE
matched "or", DOMPERIDONE matched "DOM", etc.). This script applies a
curated whitelist of genuine PDSP ligand names for Shulgin compounds,
plus a manual mapping from PDSP names to canonical Shulgin names.

Outputs:
  analysis/shulgin_clean.csv          - cleaned dataset, one row per measurement
  analysis/shulgin_compound_list.csv  - one row per compound with receptor coverage
  analysis/02_clean_report.txt        - plain-text summary
"""

import pandas as pd
import numpy as np

log = open("analysis/02_clean_report.txt", "w")
def p(*args, **kwargs):
    print(*args, **kwargs, file=log)

# ── 1. Load hits from script 01 ───────────────────────────────────────────────

df = pd.read_csv("analysis/shulgin_hits.csv", low_memory=False)
p(f"Loaded {len(df)} records, {df['ligandname'].nunique()} unique names")

# ── 2. Curated whitelist: PDSP name -> canonical Shulgin name ─────────────────
# Only include compounds that are genuinely from PIHKAL or TIHKAL,
# or are closely related reference compounds (serotonin, tryptamine).
# Excludes: antipsychotics (risperidone, ziprasidone, iloperidone),
#           antidepressants (trazodone), other drugs matched by substring.

NAME_MAP = {
    # Serotonin and tryptamine (reference compounds)
    "5-Hydroxy Tryptamine":                 "serotonin",
    "SEROTONIN":                            "serotonin",
    "TRYPTAMINE":                           "tryptamine",
    "tryptamine":                           "tryptamine",

    # DMT and analogues (TIHKAL)
    "DMT":                                  "DMT",
    "DMT,5-OH":                             "bufotenine",
    "DMT,5-MeO":                            "5-MeO-DMT",
    "DMT,5-Me":                             "5-Me-DMT",
    "DMT,4-OH":                             "psilocin",
    "DMT,4-MeO":                            "4-MeO-DMT",
    "DMT,4-AcO":                            "4-AcO-DMT",
    "DMT,7-OH":                             "7-OH-DMT",
    "DMT,7-MeO":                            "7-MeO-DMT",
    "DMT,7-Br":                             "7-Br-DMT",
    "DMT,6-MeO":                            "6-MeO-DMT",
    "DMT,5-MeO-7-M":                        "5-MeO-7-Me-DMT",
    "DMT,2-Me":                             "2-Me-DMT",

    # 5-MeO tryptamines (TIHKAL)
    "5-Methoxytryptamine":                  "5-MeO-T",
    "5-MeO-DMT":                            "5-MeO-DMT",
    "5-Methoxy-N,N-diisopropyltryptamine":  "5-MeO-DiPT",
    "5-Methoxy-N-methyl-N-isopropyltryptamine": "5-MeO-MiPT",
    "5-methoxy-N,N-diallyltryptamine":      "5-MeO-DALT",
    "5-MeOT,N,N tryptamine":               "5-MeO-T",
    "5-M-N,N-DMT":                         "5-MeO-DMT",

    # DET, DPT, DiPT (TIHKAL)
    "N,N-Diethyltryptamine":                "DET",
    "DET":                                  "DET",
    "DET,5-MeO":                            "5-MeO-DET",
    "DET,4-OH":                             "4-HO-DET",
    "N,N-Dipropyltryptamine":               "DPT",
    "N,N-Diisopropyltryptamine":            "DiPT",
    "DIPT,5-OH":                            "5-HO-DiPT",

    # AMT, alpha-methyltryptamine (TIHKAL)
    "AMT":                                  "AMT",
    "AMT,(+)":                              "AMT (+)",
    "AMT,(-)":                              "AMT (-)",
    "TRYPTAMINE,a-Me":                      "AMT",
    "TRYPTAMINE,a-Me(+)":                   "AMT (+)",
    "TRYPTAMINE,a-Me(-)":                   "AMT (-)",
    "TRYPTAMINE,a-MeISO":                   "AMT-iso",
    "TRYPTAMINE,5-MeO,a-":                  "5-MeO-AMT",

    # Psilocybin / psilocin (TIHKAL)
    "PSILOCYBIN":                           "psilocybin",
    "psilocybin":                           "psilocybin",
    "4-Hydroxy-N,N-dimethyltryptamine":     "psilocin",
    "psilocin":                             "psilocin",

    # Bufotenine (TIHKAL)
    "bufotenine":                           "bufotenine",
    "BUFOTENINE":                           "bufotenine",

    # Melatonin (reference)
    "MELATONIN":                            "melatonin",

    # LSD and analogues (TIHKAL)
    "LSD":                                  "LSD",
    "LSD,(+)":                              "LSD (+)",
    "LSD,(-)":                              "LSD (-)",
    "LSD,2-Bromo":                          "2-Bromo-LSD",
    "LSD,2-iodo":                           "2-Iodo-LSD",
    "LSD, IODO":                            "Iodo-LSD",
    "[3H] LSD":                             "LSD",

    # Ibogaine (TIHKAL)
    "IBOGAINE":                             "ibogaine",
    "ibogaine":                             "ibogaine",

    # MDA / MDMA series (PIHKAL)
    "MDMA":                                 "MDMA",
    "MDA":                                  "MDA",
    "MDA,R(-)":                             "MDA R(-)",
    "MDA, (R,S)":                           "MDA (R,S)",
    "MDE":                                  "MDE",

    # DOx phenylisopropylamines (PIHKAL)
    "DOB":                                  "DOB",
    "DOI":                                  "DOI",
    "DOM":                                  "DOM",
    "DOM,iso":                              "DOM-iso",
    "DOET":                                 "DOET",
    "DOB,OH-":                              "OH-DOB",
    "DOB,BENZYLOXY-":                       "BenzO-DOB",

    # 2C-X phenethylamines (PIHKAL)
    "4-Bromo-2,5-dimethoxyphenethylamine":  "2C-B",
    "4-Ethyl-2,5-dimethoxyphenethylamine":  "2C-E",
    "4-Ethylthio-2,5-dimethoxyphenethylamine": "2C-T-2",
    "2,5-dimethoxyphenethylamin":           "2C-H",

    # Mescaline and TMA series (PIHKAL)
    "3,4,5-Trimethoxyphenethylamine":       "mescaline",
    "3,4,5-Trimethoxyphenethylamino":       "mescaline",
    "3,4,5-Trimethoxyphenethylamini":       "mescaline",
    "TMA,2,3,5-":                           "TMA-2,3,5",
    "TMA,2,4,5-":                           "TMA-2,4,5",
    "2,4,5-TMA":                            "TMA-2,4,5",

    # Miscellaneous tryptamines
    "6-Fluoro-N,N-dimethyltryptamine":      "6-F-DMT",
    "2-phenyl-N,N-diallyltryptamine":       "2-Ph-DALT",
    "N,N-diallyltryptamine":               "DALT",
    "5-bromo-N,N-diallyltryptamine":        "5-Br-DALT",
    "5-fluoro-N,N-diallyltryptamine":       "5-F-DALT",
    "5-methoxy-N,N-diallyltryptamine":      "5-MeO-DALT",
    "4-hydroxy-N,N-diallyltryptamine":      "4-HO-DALT",
    "4-acetoxy-N,N-diallyltryptamine":      "4-AcO-DALT",
    "5-methoxy-2-methyl-N,N-diallyltryptamine": "5-MeO-2-Me-DALT",
    "5-methoxy-2-fluoro-N,N-diallyltryptamine": "5-MeO-2-F-DALT",
    "7-ethyl-N,N-diallyltryptamine":        "7-Et-DALT",
    "DET,5-MeO":                            "5-MeO-DET",
}

# ── 3. Apply whitelist ────────────────────────────────────────────────────────

df["canonical_name"] = df["ligandname"].map(NAME_MAP)
df_clean = df[df["canonical_name"].notna()].copy()

p(f"\nAfter whitelist filter:")
p(f"  Records:           {len(df_clean)}")
p(f"  Unique compounds:  {df_clean['canonical_name'].nunique()}")
p(f"  Removed (noise):   {len(df) - len(df_clean)} records")

# ── 4. Clean Ki values ────────────────────────────────────────────────────────
# kival is stored as string; kinote = ">" means censored upper bound.
# For censored values we keep them but flag them.
# Convert to numeric, compute pKi = -log10(Ki / 1e9) = 9 - log10(Ki)
# where Ki is in nM.

df_clean["ki_nM"] = pd.to_numeric(df_clean["kival"], errors="coerce")
df_clean["censored"] = df_clean["kinote"].astype(str).str.strip() == ">"
df_clean["pKi"] = 9 - np.log10(df_clean["ki_nM"])

valid = df_clean["ki_nM"].notna()
p(f"\nKi value quality:")
p(f"  Valid numeric Ki:  {valid.sum()}")
p(f"  Non-numeric/NaN:   {(~valid).sum()}")
p(f"  Censored ('>'):    {df_clean['censored'].sum()}")
p(f"  pKi range:         {df_clean['pKi'].min():.2f} to {df_clean['pKi'].max():.2f}")

# ── 5. Filter to human species where possible ─────────────────────────────────

p(f"\nSpecies breakdown:")
p(df_clean["species"].value_counts().head(10).to_string())

# Keep human; fall back to rat if no human data for a compound/receptor
human_mask = df_clean["species"].str.upper().str.contains("HUMAN", na=False)
p(f"\nHuman records: {human_mask.sum()} / {len(df_clean)}")

# ── 6. Receptor coverage per compound ────────────────────────────────────────

coverage = (df_clean[valid]
            .groupby("canonical_name")["name"]
            .apply(lambda x: sorted(x.unique().tolist()))
            .reset_index())
coverage.columns = ["canonical_name", "receptors"]
coverage["n_receptors"] = coverage["receptors"].apply(len)
coverage["has_5HT2A"]   = coverage["receptors"].apply(lambda r: "5-HT2A" in r)
coverage["has_5HT2B"]   = coverage["receptors"].apply(lambda r: "5-HT2B" in r)
coverage = coverage.sort_values("n_receptors", ascending=False)

p(f"\nReceptor coverage per compound:")
p(f"{'Compound':<25} {'n_rec':>5} {'2A':>4} {'2B':>4}  Receptors")
p("-"*80)
for _, row in coverage.iterrows():
    p(f"  {row['canonical_name']:<23} {row['n_receptors']:>5} "
      f"{'Y' if row['has_5HT2A'] else '-':>4} "
      f"{'Y' if row['has_5HT2B'] else '-':>4}  "
      f"{', '.join(row['receptors'])}")

p(f"\nCompounds with 5-HT2A data: {coverage['has_5HT2A'].sum()}")
p(f"Compounds with >= 3 receptors: {(coverage['n_receptors'] >= 3).sum()}")
p(f"Compounds with >= 5 receptors: {(coverage['n_receptors'] >= 5).sum()}")

# ── 7. Save ───────────────────────────────────────────────────────────────────

df_clean.to_csv("analysis/shulgin_clean.csv", index=False)
coverage.to_csv("analysis/shulgin_compound_list.csv", index=False)

p(f"\nSaved analysis/shulgin_clean.csv ({len(df_clean)} records)")
p(f"Saved analysis/shulgin_compound_list.csv ({len(coverage)} compounds)")

log.close()
print(f"Done. Report written to analysis/02_clean_report.txt")
print(f"Clean records: {len(df_clean)}, compounds: {df_clean['canonical_name'].nunique()}")
