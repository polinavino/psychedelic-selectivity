"""
Script 01: Inspect PDSP data and identify Shulgin compounds.
All output written to analysis/01_inspect_report.txt — nothing printed to shell.

Outputs:
  analysis/serotonin_receptors.csv
  analysis/shulgin_hits.csv
  analysis/01_inspect_report.txt
"""

import pandas as pd
import sys

# Redirect all output to report file
report_path = "analysis/01_inspect_report.txt"
log = open(report_path, "w")

def p(*args, **kwargs):
    print(*args, **kwargs, file=log)

# ── 1. Load ───────────────────────────────────────────────────────────────────

df = pd.read_csv("analysis/KiDatabase.csv", low_memory=False)
p(f"Loaded {len(df)} records")
p(f"Unique receptors: {df['name'].nunique()}")

# ── 2. Serotonin receptors ────────────────────────────────────────────────────

SEROTONIN_RECEPTORS = {
    "5-HT1A", "5-HT1B", "5-HT1D", "5-HT1E", "5-HT1F",
    "5-HT2A", "5-HT2B", "5-HT2C",
    "5-HT3", "5-HT4", "5-HT5A", "5-HT6", "5-HT7",
}

df_ht = df[df["name"].isin(SEROTONIN_RECEPTORS)].copy()
p(f"\nSerotonin receptor records: {len(df_ht)}")
p(df_ht["name"].value_counts().to_string())
df_ht.to_csv("analysis/serotonin_receptors.csv", index=False)

# ── 3. Shulgin compound list ──────────────────────────────────────────────────

PIHKAL_COMPOUNDS = {
    "2C-B", "2C-C", "2C-D", "2C-E", "2C-F", "2C-G", "2C-H",
    "2C-I", "2C-N", "2C-O", "2C-P", "2C-T", "2C-T-2", "2C-T-4",
    "2C-T-7", "2C-T-21",
    "DOB", "DOC", "DOI", "DOM", "DON", "DOET", "DOPr",
    "mescaline", "MMDA", "MMDA-2", "MMDA-3a", "MMDA-3b",
    "MDA", "MDMA", "MDE", "MDEA",
    "TMA", "TMA-2", "TMA-3", "TMA-4", "TMA-5", "TMA-6",
    "PMA", "PMMA", "4-MA",
    "escaline", "proscaline", "metaescaline",
    "Aleph", "Aleph-2", "Aleph-4", "Aleph-6", "Aleph-7",
    "ARIADNE", "BOB", "brolamfetamine",
    "G-3", "G-4", "G-5", "ganesha",
    "HOT-2", "HOT-7", "HOT-17",
    "LOPHOPHINE", "MADAM-6", "ORTHO-DOT",
}

TIHKAL_COMPOUNDS = {
    "DMT", "5-MeO-DMT", "5-MeO-DIPT", "5-MeO-DPT", "5-MeO-AMT",
    "5-MeO-MiPT", "5-MeO-DALT",
    "DET", "DPT", "DiPT", "MiPT", "DIPT",
    "AMT", "alpha-MT",
    "psilocin", "psilocybin", "baeocystin",
    "bufotenine", "bufotenin",
    "harmaline", "harmine", "tetrahydroharmine", "harmane",
    "LSD", "ALD-52", "ETH-LAD", "LSA", "LSM-775", "MLD-41", "PRO-LAD",
    "ibogaine", "noribogaine", "tabernanthine",
    "4-HO-DMT", "4-HO-DET", "4-HO-DPT", "4-HO-DiPT", "4-HO-MiPT",
    "4-HO-MET", "4-HO-pyr-T",
    "4-AcO-DMT", "4-AcO-DET", "4-AcO-DIPT", "4-AcO-MiPT",
    "5-HO-DMT", "5-OH-DMT",
    "tryptamine", "serotonin", "melatonin",
    "NET", "NMT", "NBT",
}

ALL_SHULGIN = PIHKAL_COMPOUNDS | TIHKAL_COMPOUNDS
p(f"\nShulgin compounds searched: {len(ALL_SHULGIN)}")

# ── 4. Search PDSP ────────────────────────────────────────────────────────────

df_ht["ligandname_lower"] = df_ht["ligandname"].str.lower().str.strip()
shulgin_lower = {s.lower(): s for s in ALL_SHULGIN}

exact_mask = df_ht["ligandname_lower"].isin(shulgin_lower.keys())
substr_mask = df_ht["ligandname_lower"].apply(
    lambda n: any(s in n for s in shulgin_lower.keys()) if isinstance(n, str) else False
)
df_substr_only = df_ht[substr_mask & ~exact_mask].copy()
df_hits = pd.concat([df_ht[exact_mask], df_substr_only], ignore_index=True)

p(f"\nExact matches:     {exact_mask.sum()} records")
p(f"Substring matches: {(substr_mask & ~exact_mask).sum()} records")
p(f"Total records:     {len(df_hits)}")

found_compounds = df_hits["ligandname"].value_counts()
p(f"\nUnique Shulgin compounds found in PDSP: {len(found_compounds)}")
p(found_compounds.to_string())

# ── 5. 5-HT2A coverage ───────────────────────────────────────────────────────

df_2a = df_hits[df_hits["name"] == "5-HT2A"]
found_2a = df_2a["ligandname"].value_counts()
p(f"\nWith 5-HT2A data: {len(found_2a)} compounds")
p(found_2a.to_string())

# ── 6. Receptor coverage per compound ────────────────────────────────────────

coverage = (df_hits.groupby("ligandname")["name"]
            .apply(lambda x: sorted(x.unique().tolist()))
            .reset_index())
coverage.columns = ["ligandname", "receptors"]
coverage["n_receptors"] = coverage["receptors"].apply(len)
coverage = coverage.sort_values("n_receptors", ascending=False)

p(f"\nReceptor coverage per Shulgin compound (all):")
for _, row in coverage.iterrows():
    p(f"  {row['ligandname']:<30} ({row['n_receptors']:>2} receptors): "
      f"{', '.join(row['receptors'])}")

# ── 7. Ki value check ────────────────────────────────────────────────────────

p(f"\nKi value sample (first 20):")
p(df_hits[["ligandname", "name", "kival", "kinote"]].head(20).to_string(index=False))

p(f"\nKi note values (what modifiers exist):")
p(df_hits["kinote"].value_counts().to_string())

# ── 8. Save ───────────────────────────────────────────────────────────────────

df_hits.to_csv("analysis/shulgin_hits.csv", index=False)
p(f"\nSaved analysis/shulgin_hits.csv ({len(df_hits)} records)")

log.close()
print(f"Done. Report written to {report_path}")
print(f"Found {len(found_compounds)} Shulgin compounds, "
      f"{len(found_2a)} with 5-HT2A data.")
