"""
Script 02: Prepare Shulgin compound ligands for AutoDock Vina docking.

Takes SMILES from the manual dictionary (validated in psychedelic-selectivity
analysis) for each compound with 5-HT2A pKi data in PDSP, generates 3D
conformers using OpenBabel, and saves as PDBQT files for Vina.

OpenBabel pipeline per ligand:
  SMILES -> 3D geometry (--gen3d, MMFF94 force field) -> PDBQT
  -p 7.4  : protonate at physiological pH
  --gen3d : generate 3D coordinates from scratch
  -h      : add hydrogens

Compounds without valid SMILES are skipped and logged.

Outputs:
  data/ligands/<compound_name>.pdbqt   - one per compound
  data/ligands/ligand_manifest.csv     - compound, SMILES, status
  analysis/02_ligand_prep_report.txt

Usage:
  cd psychedelic-selectivity/docking/
  python analysis/02_prepare_ligands.py
"""

import subprocess
import os
import re
import pandas as pd

os.makedirs("data/ligands", exist_ok=True)

log = open("analysis/02_ligand_prep_report.txt", "w")
def p(*args, **kwargs):
    print(*args, **kwargs, file=log)

# ── SMILES dictionary ─────────────────────────────────────────────────────────
# Validated SMILES from the psychedelic-selectivity/shulgin-selectivity
# analysis (script 04_bias_predictions.py). All checked with RDKit.

SMILES = {
    # Reference compounds
    "serotonin":        "NCCc1c[nH]c2cccc(O)c12",
    "tryptamine":       "NCCc1c[nH]c2ccccc12",
    "melatonin":        "COc1ccc2[nH]cc(CCNC(C)=O)c2c1",
    # DMT and simple analogues
    "DMT":              "CN(C)CCc1c[nH]c2ccccc12",
    "5-MeO-DMT":        "CN(C)CCc1c[nH]c2cc(OC)ccc12",
    "bufotenine":       "CN(C)CCc1c[nH]c2cc(O)ccc12",
    "psilocin":         "CN(C)CCc1[nH]c2cccc(O)c2c1",
    "5-Me-DMT":         "CN(C)CCc1c[nH]c2cc(C)ccc12",
    "6-F-DMT":          "CN(C)CCc1c[nH]c2ccc(F)cc12",
    # LSD and analogues
    "LSD":              "CCN(CC)C(=O)[C@H]1CN(C)C[C@@H]2Cc3c[nH]c4cccc(c34)[C@@H]12",
    "LSD (+)":          "CCN(CC)C(=O)[C@H]1CN(C)C[C@@H]2Cc3c[nH]c4cccc(c34)[C@@H]12",
    "2-Bromo-LSD":      "CCN(CC)C(=O)[C@H]1CN(C)C[C@@H]2Cc3c(Br)[nH]c4cccc(c34)[C@@H]12",
    # Psilocybin
    "psilocybin":       "CN(C)CCc1[nH]c2cccc(OP(=O)(O)O)c2c1",
    # Other tryptamines
    "DPT":              "CCCN(CCC)CCc1c[nH]c2ccccc12",
    "DiPT":             "CC(C)N(CC(C)C)CCc1c[nH]c2ccccc12",
    "5-MeO-DiPT":       "CC(C)N(CC(C)C)CCc1c[nH]c2cc(OC)ccc12",
    "5-MeO-MiPT":       "CN(CC(C)C)CCc1c[nH]c2cc(OC)ccc12",
    "5-MeO-T":          "NCCc1c[nH]c2cc(OC)ccc12",
    "AMT":              "CC(N)Cc1c[nH]c2ccccc12",
    "AMT (+)":          "[C@@H](Cc1c[nH]c2ccccc12)(N)C",
    "AMT (-)":          "[C@H](Cc1c[nH]c2ccccc12)(N)C",
    # DALT series
    "DALT":             "C=CCN(CC=C)CCc1c[nH]c2ccccc12",
    "5-MeO-DALT":       "C=CCN(CC=C)CCc1c[nH]c2cc(OC)ccc12",
    "5-F-DALT":         "C=CCN(CC=C)CCc1c[nH]c2cc(F)ccc12",
    "5-Br-DALT":        "C=CCN(CC=C)CCc1c[nH]c2cc(Br)ccc12",
    "4-HO-DALT":        "C=CCN(CC=C)CCc1[nH]c2cccc(O)c2c1",
    "4-AcO-DALT":       "C=CCN(CC=C)CCc1[nH]c2cccc(OC(C)=O)c2c1",
    "2-Ph-DALT":        "C=CCN(CC=C)CCc1c(-c2ccccc2)[nH]c2ccccc12",
    "7-Et-DALT":        "C=CCN(CC=C)CCc1c[nH]c2cccc(CC)c12",
    "5-MeO-2-Me-DALT":  "C=CCN(CC=C)CCc1[nH]c2cc(OC)ccc2c1C",
    "5-MeO-2-F-DALT":   "C=CCN(CC=C)CCc1[nH]c2cc(OC)ccc2c1F",
    # Ibogaine
    "ibogaine":         "COc1ccc2[nH]c3c(c2c1)C[C@H]1CCN2CC[C@@H]([C@H]12)CC3",
    # Mescaline
    "mescaline":        "COc1cc(CCN)cc(OC)c1OC",
    # MDMA/MDA
    "MDMA":             "CNC(C)Cc1ccc2c(c1)OCO2",
    "MDA":              "CC(N)Cc1ccc2c(c1)OCO2",
    "MDA R(-)":         "[C@@H](Cc1ccc2c(c1)OCO2)(N)C",
    "MDA (R,S)":        "CC(N)Cc1ccc2c(c1)OCO2",
    # DOx phenylisopropylamines
    "DOB":              "CC(N)Cc1cc(OC)c(Br)cc1OC",
    "DOI":              "CC(N)Cc1cc(OC)c(I)cc1OC",
    "DOM":              "CC(N)Cc1cc(OC)c(C)cc1OC",
    "DOET":             "CC(N)Cc1cc(OC)c(CC)cc1OC",
    # 2C-X phenethylamines
    "2C-B":             "NCCc1cc(OC)c(Br)cc1OC",
    "2C-E":             "NCCc1cc(OC)c(CC)cc1OC",
    "2C-T-2":           "NCCc1cc(OC)c(SCC)cc1OC",
    "2C-H":             "NCCc1cc(OC)ccc1OC",
}

def safe_name(compound):
    """Convert compound name to filesystem-safe string."""
    return re.sub(r'[^\w\-]', '_', compound).strip('_')

# ── Prepare each ligand ───────────────────────────────────────────────────────

p("Ligand Preparation Report")
p("=" * 60)
p(f"Total compounds: {len(SMILES)}")
p("")

manifest_rows = []
success, failed = 0, 0

for compound, smiles in SMILES.items():
    fname   = safe_name(compound)
    out_pdbqt = f"data/ligands/{fname}.pdbqt"

    cmd = [
        "obabel",
        f"-:{smiles}",     # SMILES input via stdin
        "-O", out_pdbqt,
        "--gen3d",          # generate 3D coordinates
        "-p", "7.4",        # protonate at pH 7.4
        "-h",               # add hydrogens
        "--ff", "MMFF94",   # MMFF94 force field for minimisation
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if os.path.exists(out_pdbqt) and os.path.getsize(out_pdbqt) > 0:
        n_atoms = sum(1 for l in open(out_pdbqt)
                      if l.startswith(("ATOM", "HETATM")))
        status  = "OK"
        success += 1
        p(f"  OK  {compound:<28} {n_atoms:>3} atoms -> {fname}.pdbqt")
    else:
        status = "FAILED"
        failed += 1
        p(f"  FAIL {compound:<28} {result.stderr.strip()[:60]}")

    manifest_rows.append({
        "compound":  compound,
        "safe_name": fname,
        "smiles":    smiles,
        "pdbqt":     out_pdbqt if status == "OK" else "",
        "status":    status,
    })

manifest = pd.DataFrame(manifest_rows)
manifest.to_csv("data/ligands/ligand_manifest.csv", index=False)

p(f"\nSummary: {success} succeeded, {failed} failed")
p(f"Manifest saved: data/ligands/ligand_manifest.csv")

log.close()
print(f"Done. {success}/{len(SMILES)} ligands prepared.")
print(f"Report: analysis/02_ligand_prep_report.txt")
print(f"Failed: {failed}")
