"""
Script 01: Explore ChEMBL functional assay data for 5-HT2A biased agonism.

Compounds are only considered paired if they appear in both a Gq-pathway
assay and a beta-arrestin assay from the SAME paper (document_chembl_id),
ensuring assay conditions are comparable.

Outputs:
  analysis/assay_summary.csv       - all functional assays with pathway labels
  analysis/gq_assays.csv           - Gq-classified assays
  analysis/barr_assays.csv         - beta-arrestin-classified assays
  analysis/paired_documents.csv    - papers containing both pathway types
  analysis/paired_compounds.csv    - compounds paired within the same paper
  analysis/01_explore_report.txt   - plain-text summary
"""

import requests
import pandas as pd
import re

CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"
TARGET_CHEMBL_ID = "CHEMBL224"

def fetch_all_pages(url, params, key):
    records = []
    params = dict(params)
    params["limit"] = 1000
    params["offset"] = 0
    while True:
        r = requests.get(url, params=params)
        r.raise_for_status()
        data = r.json()
        records.extend(data.get(key, []))
        if data.get("page_meta", {}).get("next") is None:
            break
        params["offset"] += params["limit"]
        print(f"  fetched {len(records)} so far...")
    return records

print("Fetching functional assay activity records for CHEMBL224...")
activity_records = fetch_all_pages(
    f"{CHEMBL_API}/activity",
    {"target_chembl_id": TARGET_CHEMBL_ID, "assay_type": "F", "format": "json"},
    "activities"
)
print(f"Total functional activity records: {len(activity_records)}")

df = pd.DataFrame(activity_records)

keep = [c for c in ["assay_chembl_id", "assay_description", "document_chembl_id",
                     "standard_type", "standard_units", "molecule_chembl_id",
                     "standard_value", "standard_relation"] if c in df.columns]
df = df[keep].copy()
df["standard_type"] = df["standard_type"].str.strip().str.upper()

# ── Assay-level summary ───────────────────────────────────────────────────────

assay_info = (
    df.groupby(["assay_chembl_id", "assay_description", "document_chembl_id"])
    .size().reset_index(name="n_records")
    .sort_values("n_records", ascending=False)
)

GQ_KEYWORDS = [
    r"gq", r"g[- ]?protein", r"g[- ]?alpha[- ]?q",
    r"ip[1-3]", r"inositol", r"calcium", r"ca2\+",
    r"phospholipase", r"ip accumulation", r"cre[- ]?luc",
    r"aequorin", r"fluo[- ]?4", r"fluo4", r"miniG", r"mini[- ]?g",
]
BARR_KEYWORDS = [
    r"beta[- ]?arrestin", r"b[- ]?arrestin", r"\barr2?\b",
    r"barrestin", r"bret.*arrestin", r"arrestin.*bret",
    r"tango", r"path[- ]?hunter",
]

def classify_assay(desc):
    if not isinstance(desc, str):
        return "unknown"
    d = desc.lower()
    gq   = any(re.search(p, d) for p in GQ_KEYWORDS)
    barr = any(re.search(p, d) for p in BARR_KEYWORDS)
    if barr and gq:  return "both"
    if barr:         return "beta-arrestin"
    if gq:           return "gq"
    return "unknown"

assay_info["pathway"] = assay_info["assay_description"].apply(classify_assay)
assay_info.to_csv("analysis/assay_summary.csv", index=False)

gq_assays   = assay_info[assay_info["pathway"] == "gq"]
barr_assays = assay_info[assay_info["pathway"] == "beta-arrestin"]
both_assays = assay_info[assay_info["pathway"] == "both"]
gq_assays.to_csv("analysis/gq_assays.csv", index=False)
barr_assays.to_csv("analysis/barr_assays.csv", index=False)

print(f"\nUnique assays: {len(assay_info)}")
print(assay_info["pathway"].value_counts().to_string())

# ── Document-level pairing ────────────────────────────────────────────────────

gq_docs   = set(gq_assays["document_chembl_id"]) | set(both_assays["document_chembl_id"])
barr_docs = set(barr_assays["document_chembl_id"]) | set(both_assays["document_chembl_id"])
paired_docs = gq_docs & barr_docs
print(f"\nDocuments with Gq assays:          {len(gq_docs)}")
print(f"Documents with barr assays:        {len(barr_docs)}")
print(f"Documents with BOTH pathway types: {len(paired_docs)}")

paired_doc_rows = assay_info[assay_info["document_chembl_id"].isin(paired_docs)]
paired_doc_rows.to_csv("analysis/paired_documents.csv", index=False)

# ── Compound pairing within documents ────────────────────────────────────────

gq_assay_ids   = set(gq_assays["assay_chembl_id"]) | set(both_assays["assay_chembl_id"])
barr_assay_ids = set(barr_assays["assay_chembl_id"]) | set(both_assays["assay_chembl_id"])

emax_types = {"EMAX", "EFFICACY", "INTRINSIC ACTIVITY", "% ACTIVITY",
              "% EFFICACY", "% EFFECT", "RELATIVE EFFICACY"}
ec50_types = {"EC50", "PEC50", "LOG EC50"}

df_paired_docs = df[df["document_chembl_id"].isin(paired_docs)].copy()

df_gq   = df_paired_docs[df_paired_docs["assay_chembl_id"].isin(gq_assay_ids)]
df_barr = df_paired_docs[df_paired_docs["assay_chembl_id"].isin(barr_assay_ids)]

compounds_gq_emax   = set(df_gq[df_gq["standard_type"].isin(emax_types)]["molecule_chembl_id"].dropna())
compounds_barr_emax = set(df_barr[df_barr["standard_type"].isin(emax_types)]["molecule_chembl_id"].dropna())
paired_emax = compounds_gq_emax & compounds_barr_emax

compounds_gq_ec50   = set(df_gq[df_gq["standard_type"].isin(ec50_types)]["molecule_chembl_id"].dropna())
compounds_barr_ec50 = set(df_barr[df_barr["standard_type"].isin(ec50_types)]["molecule_chembl_id"].dropna())
paired_ec50 = compounds_gq_ec50 & compounds_barr_ec50

both_paired = paired_emax & paired_ec50

print(f"\nWithin paired documents:")
print(f"  Compounds with Emax in both pathways:    {len(paired_emax)}")
print(f"  Compounds with EC50 in both pathways:    {len(paired_ec50)}")
print(f"  Compounds with Emax AND EC50 in both:    {len(both_paired)}")

print("\nStandard types in Gq assays (paired docs):")
print(df_gq["standard_type"].value_counts().head(10).to_string())
print("\nStandard types in barr assays (paired docs):")
print(df_barr["standard_type"].value_counts().head(10).to_string())

if paired_emax:
    rows_gq   = df_gq[df_gq["standard_type"].isin(emax_types) &
                       df_gq["molecule_chembl_id"].isin(paired_emax)].copy()
    rows_barr = df_barr[df_barr["standard_type"].isin(emax_types) &
                         df_barr["molecule_chembl_id"].isin(paired_emax)].copy()
    rows_gq["pathway"]   = "gq"
    rows_barr["pathway"] = "beta-arrestin"
    pd.concat([rows_gq, rows_barr], ignore_index=True).to_csv(
        "analysis/paired_compounds.csv", index=False)
    print(f"\nSaved analysis/paired_compounds.csv")

print("\nSample Gq assay descriptions (paired docs):")
for desc in gq_assays[gq_assays["document_chembl_id"].isin(paired_docs)]["assay_description"].head(5):
    print(f"  - {desc}")
print("\nSample barr assay descriptions (paired docs):")
for desc in barr_assays[barr_assays["document_chembl_id"].isin(paired_docs)]["assay_description"].head(5):
    print(f"  - {desc}")

report = f"""
5-HT2A Biased Agonism Data Exploration - ChEMBL (CHEMBL224)
============================================================
Compounds are paired only within the same document (paper).

Total functional activity records:   {len(activity_records)}
Unique functional assays:            {len(assay_info)}

Assay pathway classification:
  Gq-pathway:       {len(gq_assays)}
  Beta-arrestin:    {len(barr_assays)}
  Both keywords:    {len(both_assays)}
  Unclassified:     {(assay_info['pathway'] == 'unknown').sum()}

Documents:
  With Gq assays:          {len(gq_docs)}
  With barr assays:        {len(barr_docs)}
  With BOTH (qualifying):  {len(paired_docs)}

Within qualifying documents:
  Compounds with Emax in both pathways:      {len(paired_emax)}
  Compounds with EC50 in both pathways:      {len(paired_ec50)}
  Compounds with Emax AND EC50 in both:      {len(both_paired)}

Verdict:
  {"Good - proceed to bias factor computation." if len(paired_emax) >= 30 else "Small dataset (<30). ChEMBL alone may be insufficient; consider literature curation."}

Gq assay descriptions (sample, paired docs):
{chr(10).join("  - " + str(d) for d in gq_assays[gq_assays["document_chembl_id"].isin(paired_docs)]["assay_description"].head(8))}

Barr assay descriptions (sample, paired docs):
{chr(10).join("  - " + str(d) for d in barr_assays[barr_assays["document_chembl_id"].isin(paired_docs)]["assay_description"].head(8))}
"""
print(report)
with open("analysis/01_explore_report.txt", "w") as f:
    f.write(report)
print("Done.")
