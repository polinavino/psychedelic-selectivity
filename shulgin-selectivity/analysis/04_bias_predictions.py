"""
Script 04: Predict signalling bias at 5-HT2A for Shulgin compounds.

Applies the Lasso QSAR model trained on the biased agonism dataset
(psychedelic-selectivity/analysis/combined_dataset.csv) to predict
log_bias = log10((Emax_barr/EC50_barr)/(Emax_Gq/EC50_Gq)) for each
Shulgin compound with a SMILES string.

log_bias > 0  : predicted beta-arrestin biased (potentially non-hallucinogenic)
log_bias < 0  : predicted Gq biased (potentially hallucinogenic)
log_bias ~ 0  : predicted balanced

IMPORTANT CAVEAT: The QSAR model showed R²=0.78 in LOO-CV but R²<0 in
leave-one-series-out CV, indicating it learned scaffold identity rather
than generalizable structure-activity relationships. Predictions for
Shulgin compounds (a new chemical series) should be treated as exploratory
hypotheses only, not reliable quantitative predictions.

Parameters (change here to adjust behaviour):
  RADIUS        : Morgan fingerprint radius (2 = ECFP4)
  N_BITS        : fingerprint length in bits
  LASSO_ALPHA   : Lasso regularisation strength
  BIAS_THRESHOLD: |log_bias| above which a compound is classified as biased

Outputs:
  analysis/bias_predictions.csv
  analysis/04_bias_report.txt
  analysis/bias_vs_selectivity.png
"""

import pandas as pd
import numpy as np
import requests
import time
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from rdkit import Chem
from rdkit.Chem import AllChem

from sklearn.linear_model import Lasso
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

# ── Parameters ────────────────────────────────────────────────────────────────

RADIUS         = 2      # Morgan fingerprint radius (ECFP4)
N_BITS         = 2048   # fingerprint bit length
LASSO_ALPHA    = 0.01   # Lasso regularisation (matches training in script 05)
BIAS_THRESHOLD = 0.3    # |log_bias| threshold for barr/Gq classification

TRAINING_DATA  = "../analysis/combined_dataset.csv"

# ── Setup ─────────────────────────────────────────────────────────────────────

log = open("analysis/04_bias_report.txt", "w")
def p(*args, **kwargs):
    print(*args, **kwargs, file=log)

# ── 1. Load training data ─────────────────────────────────────────────────────

try:
    train_df = pd.read_csv(TRAINING_DATA)
    p(f"Loaded training data: {len(train_df)} compounds from {TRAINING_DATA}")
except FileNotFoundError:
    p(f"ERROR: Could not find {TRAINING_DATA}")
    p("Run from shulgin-selectivity/ directory, with the main repo at ../")
    raise

train_df = train_df[train_df["log_bias"].notna() & train_df["smiles"].notna()]
p(f"Training compounds with log_bias and SMILES: {len(train_df)}")

# ── 2. Training fingerprints ──────────────────────────────────────────────────

def get_fp(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return np.array(AllChem.GetMorganFingerprintAsBitVect(
        mol, radius=RADIUS, nBits=N_BITS))

train_fps, train_valid = [], []
for i, row in train_df.iterrows():
    fp = get_fp(row["smiles"])
    if fp is not None:
        train_fps.append(fp)
        train_valid.append(i)

X_train    = np.array(train_fps)
train_df   = train_df.loc[train_valid].reset_index(drop=True)
y_train    = train_df["log_bias"].values
nonzero    = X_train.var(axis=0) > 0
X_train_nz = X_train[:, nonzero]
p(f"Training fingerprints: {X_train_nz.shape} ({nonzero.sum()} informative bits)")

# ── 3. Train Lasso on full training set ───────────────────────────────────────

model = Pipeline([
    ("scaler", StandardScaler()),
    ("lasso",  Lasso(alpha=LASSO_ALPHA, max_iter=5000))
])
model.fit(X_train_nz, y_train)
p(f"Lasso model trained (alpha={LASSO_ALPHA})")
p(f"Training log_bias range: {y_train.min():.3f} to {y_train.max():.3f}")

# ── 4. SMILES for Shulgin compounds ───────────────────────────────────────────
# All manually curated and validated. PubChem lookup is attempted first
# for any compound not in this dictionary; the dictionary takes priority.
#
# Key corrections vs naive naming:
#   psilocin    : 4-OH-DMT (C4 of indole), not 5-OH (bufotenine)
#   bufotenine  : 5-OH-DMT (C5 of indole)
#   MDA         : alpha-methylphenethylamine scaffold (amphetamine), not phenethylamine
#   DiPT        : N,N-diisopropyltryptamine — both N-substituents are isopropyl
#   5-MeO-DiPT  : same as DiPT but with 5-methoxy
#   5-MeO-MiPT  : N-methyl-N-isopropyl (one methyl, one isopropyl)
#   2C-H        : 2,5-dimethoxyphenethylamine (no 4-substituent)
#   4-HO-DALT   : 4-hydroxy on indole C4 (adjacent to C3 where ethylamine attaches)
#   4-AcO-DALT  : 4-acetoxy on indole C4
#   2-Ph-DALT   : 2-phenyltryptamine N,N-diallyl (phenyl at C2 of indole)

MANUAL_SMILES = {
    # Reference compounds
    "serotonin":        "NCCc1c[nH]c2cccc(O)c12",
    "tryptamine":       "NCCc1c[nH]c2ccccc12",
    "melatonin":        "COc1ccc2[nH]cc(CCNC(C)=O)c2c1",

    # DMT and simple analogues (TIHKAL)
    "DMT":              "CN(C)CCc1c[nH]c2ccccc12",
    "5-MeO-DMT":        "CN(C)CCc1c[nH]c2cc(OC)ccc12",
    "bufotenine":       "CN(C)CCc1c[nH]c2cc(O)ccc12",       # 5-OH-DMT
    "psilocin":         "CN(C)CCc1[nH]c2cccc(O)c2c1",       # 4-OH-DMT
    "5-Me-DMT":         "CN(C)CCc1c[nH]c2cc(C)ccc12",
    "6-F-DMT":          "CN(C)CCc1c[nH]c2ccc(F)cc12",

    # LSD and analogues (TIHKAL)
    "LSD":              "CCN(CC)C(=O)[C@H]1CN(C)C[C@@H]2Cc3c[nH]c4cccc(c34)[C@@H]12",
    "LSD (+)":          "CCN(CC)C(=O)[C@H]1CN(C)C[C@@H]2Cc3c[nH]c4cccc(c34)[C@@H]12",
    "2-Bromo-LSD":      "CCN(CC)C(=O)[C@H]1CN(C)C[C@@H]2Cc3c(Br)[nH]c4cccc(c34)[C@@H]12",

    # Psilocybin (TIHKAL)
    "psilocybin":       "CN(C)CCc1[nH]c2cccc(OP(=O)(O)O)c2c1",

    # Other tryptamines (TIHKAL)
    "DPT":              "CCCN(CCC)CCc1c[nH]c2ccccc12",
    "DiPT":             "CC(C)N(CC(C)C)CCc1c[nH]c2ccccc12",
    "5-MeO-DiPT":       "CC(C)N(CC(C)C)CCc1c[nH]c2cc(OC)ccc12",
    "5-MeO-MiPT":       "CN(CC(C)C)CCc1c[nH]c2cc(OC)ccc12",
    "5-MeO-T":          "NCCc1c[nH]c2cc(OC)ccc12",
    "AMT":              "CC(N)Cc1c[nH]c2ccccc12",
    "AMT (+)":          "[C@@H](Cc1c[nH]c2ccccc12)(N)C",
    "AMT (-)":          "[C@H](Cc1c[nH]c2ccccc12)(N)C",

    # DALT series (TIHKAL)
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

    # Ibogaine (TIHKAL)
    "ibogaine":         "COc1ccc2[nH]c3c(c2c1)C[C@H]1CCN2CC[C@@H]([C@H]12)CC3",

    # Mescaline (PIHKAL)
    "mescaline":        "COc1cc(CCN)cc(OC)c1OC",

    # MDMA / MDA series (PIHKAL)
    "MDMA":             "CNC(C)Cc1ccc2c(c1)OCO2",
    "MDA":              "CC(N)Cc1ccc2c(c1)OCO2",
    "MDA R(-)":         "[C@@H](Cc1ccc2c(c1)OCO2)(N)C",
    "MDA (R,S)":        "CC(N)Cc1ccc2c(c1)OCO2",

    # DOx phenylisopropylamines (PIHKAL)
    "DOB":              "CC(N)Cc1cc(OC)c(Br)cc1OC",
    "DOI":              "CC(N)Cc1cc(OC)c(I)cc1OC",
    "DOM":              "CC(N)Cc1cc(OC)c(C)cc1OC",
    "DOET":             "CC(N)Cc1cc(OC)c(CC)cc1OC",

    # 2C-X phenethylamines (PIHKAL)
    "2C-B":             "NCCc1cc(OC)c(Br)cc1OC",
    "2C-E":             "NCCc1cc(OC)c(CC)cc1OC",
    "2C-T-2":           "NCCc1cc(OC)c(SCC)cc1OC",
    "2C-H":             "NCCc1cc(OC)ccc1OC",
}

# ── 5. Validate SMILES with RDKit ─────────────────────────────────────────────

p("\nSMILES validation:")
invalid = []
for name, smi in MANUAL_SMILES.items():
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        invalid.append(name)
        p(f"  INVALID: {name}: {smi}")
    else:
        p(f"  OK: {name}")

if invalid:
    p(f"\nWARNING: {len(invalid)} invalid SMILES — these compounds will be skipped")
else:
    p(f"\nAll {len(MANUAL_SMILES)} SMILES valid.")

# ── 6. Load Shulgin compounds and assign SMILES ───────────────────────────────

selectivity = pd.read_csv("analysis/selectivity_profiles.csv")
compounds   = selectivity["compound"].tolist()
p(f"\nShulgin compounds to predict: {len(compounds)}")

PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

def fetch_smiles_pubchem(name):
    try:
        url = f"{PUBCHEM}/compound/name/{requests.utils.quote(name)}/property/IsomericSMILES/JSON"
        r = requests.get(url, timeout=10)
        if r.status_code == 200:
            return r.json()["PropertyTable"]["Properties"][0]["IsomericSMILES"]
    except Exception:
        pass
    return None

smiles_dict = {}
p("\nSMILES source per compound:")
for name in compounds:
    if name in MANUAL_SMILES:
        smi = MANUAL_SMILES[name]
        if Chem.MolFromSmiles(smi) is not None:
            smiles_dict[name] = smi
            p(f"  {name:<25} [manual]")
        else:
            smiles_dict[name] = None
            p(f"  {name:<25} [manual — INVALID, skipping]")
    else:
        smi = fetch_smiles_pubchem(name)
        if smi and Chem.MolFromSmiles(smi) is not None:
            smiles_dict[name] = smi
            p(f"  {name:<25} [PubChem] {smi[:50]}")
        else:
            smiles_dict[name] = None
            p(f"  {name:<25} NOT FOUND")
        time.sleep(0.3)

found = sum(1 for v in smiles_dict.values() if v)
p(f"\nSMILES available: {found} / {len(compounds)}")
p(f"Missing: {[n for n, s in smiles_dict.items() if s is None]}")

# ── 7. Predict bias ───────────────────────────────────────────────────────────

predictions = []
for name, smi in smiles_dict.items():
    if smi is None:
        predictions.append({"compound": name, "smiles": None,
                            "predicted_log_bias": np.nan,
                            "valid_smiles": False})
        continue
    fp = get_fp(smi)
    fp_nz = fp[nonzero].reshape(1, -1)
    pred  = model.predict(fp_nz)[0]
    predictions.append({"compound": name, "smiles": smi,
                        "predicted_log_bias": pred,
                        "valid_smiles": True})

pred_df = pd.DataFrame(predictions)

def classify(x):
    if pd.isna(x): return "no prediction"
    if x >  BIAS_THRESHOLD: return "barr-biased"
    if x < -BIAS_THRESHOLD: return "Gq-biased"
    return "balanced"

pred_df["bias_class"] = pred_df["predicted_log_bias"].apply(classify)

# Merge with selectivity profiles
result = selectivity.merge(pred_df, on="compound", how="left")
result.to_csv("analysis/bias_predictions.csv", index=False)

# ── 8. Report ─────────────────────────────────────────────────────────────────

p(f"\n{'='*72}")
p(f"PREDICTED SIGNALLING BIAS — Shulgin compounds at 5-HT2A")
p(f"{'='*72}")
p(f"\nParameters: RADIUS={RADIUS}, N_BITS={N_BITS}, LASSO_ALPHA={LASSO_ALPHA}")
p(f"Bias threshold: |log_bias| > {BIAS_THRESHOLD}")
p(f"\n*** CAVEAT: model R²=0.78 LOO-CV but R²<0 LOSO-CV — cross-series")
p(f"    predictions are exploratory only. ***\n")

p(f"{'Compound':<22} {'predLogBias':>12} {'Class':<14} {'2A-2B':>7}  Top rec")
p("-"*72)
for _, row in result.sort_values("predicted_log_bias", ascending=False).iterrows():
    lb  = row["predicted_log_bias"]
    sel = row.get("selectivity_2A_over_2B", np.nan)
    lb_str  = f"{lb:+.3f}" if pd.notna(lb) else "   n/a"
    sel_str = f"{sel:+.2f}" if pd.notna(sel) else "  n/a"
    p(f"  {row['compound']:<20} {lb_str:>12} {row['bias_class']:<14} "
      f"{sel_str:>7}  {row.get('top_receptor','?')}")

p(f"\nBias class distribution:")
p(pred_df["bias_class"].value_counts().to_string())

p(f"\nPredicted log_bias range: "
  f"{pred_df['predicted_log_bias'].min():.3f} to "
  f"{pred_df['predicted_log_bias'].max():.3f}")

# ── 9. Plot ───────────────────────────────────────────────────────────────────

colors = {"barr-biased":   "#55A868",
          "Gq-biased":     "#DD8452",
          "balanced":      "#4C72B0",
          "no prediction": "gray"}

fig, ax = plt.subplots(figsize=(10, 7))
for _, row in result.iterrows():
    x = row.get("selectivity_2A_over_2B", np.nan)
    y = row["predicted_log_bias"]
    if pd.isna(x) or pd.isna(y):
        continue
    ax.scatter(x, y, color=colors.get(row["bias_class"], "gray"),
               s=60, alpha=0.85, zorder=3)
    ax.annotate(row["compound"][:14], (x, y),
                fontsize=6, alpha=0.8,
                xytext=(4, 4), textcoords="offset points")

ax.axhline(0,               color="gray",   lw=0.8, linestyle="--")
ax.axhline( BIAS_THRESHOLD, color="green",  lw=0.5, linestyle=":")
ax.axhline(-BIAS_THRESHOLD, color="orange", lw=0.5, linestyle=":")
ax.axvline(0,               color="gray",   lw=0.8, linestyle="--")

ax.set_xlabel("5-HT2A/2B selectivity (pKi_2A - pKi_2B)\n"
              "positive = 2A-selective (safer), negative = 2B-selective (cardiac risk)")
ax.set_ylabel(f"Predicted log bias at 5-HT2A\n"
              f"positive = β-arr biased, negative = Gq biased "
              f"(threshold ±{BIAS_THRESHOLD})")
ax.set_title("Shulgin compounds: receptor selectivity vs predicted signalling bias\n"
             "(bias predictions exploratory only — LOSO R²<0)")

from matplotlib.patches import Patch
legend = [Patch(color=c, label=l) for l, c in colors.items()
          if l != "no prediction"]
ax.legend(handles=legend, fontsize=8, title="Predicted bias class")
plt.tight_layout()
plt.savefig("analysis/bias_vs_selectivity.png", dpi=150)
p("Saved analysis/bias_vs_selectivity.png")

log.close()
print("Done. Report: analysis/04_bias_report.txt")
print(f"Predictions: {pred_df['valid_smiles'].sum()} / {len(pred_df)}")
