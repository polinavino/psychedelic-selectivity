"""
Script 02: Identify source papers and compound names for the 29 paired
ChEMBL compounds (Emax + EC50 in both Gq and beta-arrestin pathways,
within the same document).

Outputs:
  analysis/source_papers.csv       - DOIs and titles of the 2 qualifying papers
  analysis/paired_compound_names.csv - compound names/SMILES for the 29 compounds
  analysis/02_sources_report.txt   - plain-text summary
"""

import requests
import pandas as pd

CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"

# ── Load paired compounds from script 01 output ──────────────────────────────

df = pd.read_csv("analysis/paired_compounds.csv")
print(f"Loaded {len(df)} paired records for {df['molecule_chembl_id'].nunique()} compounds")
print(f"From documents: {df['document_chembl_id'].unique()}")

doc_ids      = df["document_chembl_id"].dropna().unique().tolist()
compound_ids = df["molecule_chembl_id"].dropna().unique().tolist()

# ── Fetch paper metadata ──────────────────────────────────────────────────────

print("\nFetching paper metadata...")
papers = []
for doc_id in doc_ids:
    r = requests.get(f"{CHEMBL_API}/document/{doc_id}.json")
    r.raise_for_status()
    d = r.json()
    papers.append({
        "document_chembl_id": doc_id,
        "doi":                d.get("doi"),
        "title":              d.get("title"),
        "journal":            d.get("journal"),
        "year":               d.get("year"),
        "authors":            d.get("authors"),
        "pubmed_id":          d.get("pubmed_id"),
    })
    print(f"  {doc_id}: {d.get('title', 'no title')[:80]}")

df_papers = pd.DataFrame(papers)
df_papers.to_csv("analysis/source_papers.csv", index=False)

# ── Fetch compound names and SMILES ──────────────────────────────────────────

print(f"\nFetching metadata for {len(compound_ids)} compounds...")
compounds = []
for cid in compound_ids:
    r = requests.get(f"{CHEMBL_API}/molecule/{cid}.json")
    r.raise_for_status()
    d = r.json()
    props = d.get("molecule_properties") or {}
    struct = d.get("molecule_structures") or {}
    compounds.append({
        "molecule_chembl_id":  cid,
        "pref_name":           d.get("pref_name"),
        "max_phase":           d.get("max_phase"),
        "smiles":              struct.get("canonical_smiles"),
        "inchi_key":           struct.get("standard_inchi_key"),
        "mw":                  props.get("full_mw"),
        "alogp":               props.get("alogp"),
        "hbd":                 props.get("hbd"),
        "hba":                 props.get("hba"),
        "molecule_type":       d.get("molecule_type"),
    })
    name = d.get("pref_name") or cid
    print(f"  {cid}: {name}")

df_compounds = pd.DataFrame(compounds)

# Merge with paired activity data to show Emax/EC50 per pathway per compound
df_gq   = df[df["pathway"] == "gq"][["molecule_chembl_id", "document_chembl_id",
                                      "standard_type", "standard_value"]].copy()
df_barr = df[df["pathway"] == "beta-arrestin"][["molecule_chembl_id",
                                                 "standard_type", "standard_value"]].copy()

# Pivot to get Emax and EC50 per pathway as columns
def pivot_pathway(df_path, prefix):
    return (df_path.groupby(["molecule_chembl_id", "standard_type"])["standard_value"]
            .mean().unstack("standard_type")
            .rename(columns=lambda c: f"{prefix}_{c.lower()}"))

gq_pivot   = pivot_pathway(df_gq,   "gq")
barr_pivot = pivot_pathway(df_barr, "barr")

df_merged = (df_compounds
             .merge(gq_pivot,   on="molecule_chembl_id", how="left")
             .merge(barr_pivot, on="molecule_chembl_id", how="left")
             .merge(df_gq[["molecule_chembl_id", "document_chembl_id"]].drop_duplicates(),
                    on="molecule_chembl_id", how="left"))

df_merged.to_csv("analysis/paired_compound_names.csv", index=False)
print(f"\nSaved analysis/paired_compound_names.csv")

# ── Report ────────────────────────────────────────────────────────────────────

print("\n" + "="*60)
print("SOURCE PAPERS")
print("="*60)
for _, row in df_papers.iterrows():
    print(f"\n  Document:  {row['document_chembl_id']}")
    print(f"  Title:     {row['title']}")
    print(f"  Journal:   {row['journal']} ({row['year']})")
    print(f"  DOI:       {row['doi']}")
    print(f"  PubMed:    {row['pubmed_id']}")

print("\n" + "="*60)
print("COMPOUNDS")
print("="*60)
for _, row in df_merged.iterrows():
    name = row["pref_name"] or row["molecule_chembl_id"]
    gq_emax  = row.get("gq_emax",  float("nan"))
    gq_ec50  = row.get("gq_ec50",  float("nan"))
    ba_emax  = row.get("barr_emax", float("nan"))
    ba_ec50  = row.get("barr_ec50", float("nan"))
    print(f"  {name:<35} Gq Emax={gq_emax:6.1f}  Gq EC50={gq_ec50:8.2f}  "
          f"bArr Emax={ba_emax:6.1f}  bArr EC50={ba_ec50:8.2f}")

report = f"""
5-HT2A Biased Agonism - Source Papers and Compounds
====================================================

Qualifying papers (both Gq and beta-arrestin assays, BRET platform):
{chr(10).join(f"  [{row['document_chembl_id']}] {row['title']} -- {row['journal']} ({row['year']}) -- DOI: {row['doi']}" for _, row in df_papers.iterrows())}

Total paired compounds: {len(df_compounds)}
Compounds with SMILES:  {df_compounds['smiles'].notna().sum()}
Compounds with names:   {df_compounds['pref_name'].notna().sum()}

Chemical series present (based on names):
{chr(10).join("  - " + str(n) for n in df_compounds["pref_name"].dropna().tolist())}

Next steps:
  These {len(df_compounds)} compounds form the ChEMBL core dataset.
  Supplement with manual curation from the papers identified above
  and related publications to reach ~80-100 compounds total.
"""
print(report)
with open("analysis/02_sources_report.txt", "w") as f:
    f.write(report)
print("Done.")
