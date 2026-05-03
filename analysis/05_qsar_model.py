"""
Script 05: QSAR model predicting signalling bias at 5-HT2A from molecular structure.

Target: log10 bias factor (continuous regression)
Features: Morgan fingerprints (ECFP4, radius=2, 2048 bits)

Two cross-validation schemes:
  1. Leave-One-Out (LOO-CV): each compound held out in turn.
     Tests interpolation within the chemical space of the dataset.

  2. Leave-One-Series-Out (LOSO-CV): each chemical series held out
     in turn while training on the other two.
     Tests whether the model generalises ACROSS scaffolds, or whether
     it is merely learning scaffold identity (inter-series confounding).

     If LOO R² >> LOSO R², the model is likely memorising scaffold
     identity rather than learning the structure-activity relationship.

Outputs:
  analysis/05_qsar_report.txt
  analysis/qsar_predictions.png
  analysis/feature_importance.png
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

from rdkit import Chem
from rdkit.Chem import AllChem

from sklearn.linear_model import Ridge, Lasso
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.model_selection import LeaveOneOut, cross_val_predict
from sklearn.metrics import r2_score, mean_absolute_error
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

# ── 1. Load dataset ───────────────────────────────────────────────────────────

df = pd.read_csv("analysis/combined_dataset.csv")
df = df[df["log_bias"].notna() & df["smiles"].notna()].reset_index(drop=True)
print(f"Loaded {len(df)} compounds")

# ── 2. Compute Morgan fingerprints ────────────────────────────────────────────

RADIUS = 2
N_BITS = 2048

fps, valid_idx = [], []
for i, row in df.iterrows():
    mol = Chem.MolFromSmiles(row["smiles"])
    if mol is None:
        print(f"  Invalid SMILES: {row['compound_name']} — skipping")
        continue
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=RADIUS, nBits=N_BITS)
    fps.append(np.array(fp))
    valid_idx.append(i)

X_full = np.array(fps)
df = df.loc[valid_idx].reset_index(drop=True)
y = df["log_bias"].values

# Remove zero-variance bits
nonzero = X_full.var(axis=0) > 0
X = X_full[:, nonzero]
print(f"Fingerprints: {X.shape[0]} compounds x {X.shape[1]} informative bits")

# ── 3. Define models ──────────────────────────────────────────────────────────

def make_models():
    return {
        "Ridge":     Pipeline([("scaler", StandardScaler()),
                                ("model",  Ridge(alpha=1.0))]),
        "Lasso":     Pipeline([("scaler", StandardScaler()),
                                ("model",  Lasso(alpha=0.01, max_iter=5000))]),
        "RF":        RandomForestRegressor(n_estimators=200, max_features="sqrt",
                                           random_state=42, n_jobs=-1),
        "GradBoost": GradientBoostingRegressor(n_estimators=100, max_depth=3,
                                               learning_rate=0.05, random_state=42),
    }

def evaluate(y_true, y_pred, label=""):
    r2  = r2_score(y_true, y_pred)
    mae = mean_absolute_error(y_true, y_pred)
    r, p = stats.pearsonr(y_true, y_pred)
    return {"label": label, "r2": r2, "mae": mae, "r": r, "p": p,
            "y_pred": y_pred}

def print_table(results, title):
    print(f"\n{title}")
    print("="*65)
    print(f"{'Model':<12} {'R²':>8} {'MAE':>8} {'Pearson r':>10} {'p-value':>14}")
    print("="*65)
    for name, res in results.items():
        print(f"{name:<12} {res['r2']:>8.3f} {res['mae']:>8.3f} "
              f"{res['r']:>10.3f} {res['p']:>14.2e}")
    print("="*65)

# ── 4. LOO-CV ─────────────────────────────────────────────────────────────────

print("\n--- Leave-One-Out Cross-Validation ---")
loo_results = {}
for name, model in make_models().items():
    y_pred = cross_val_predict(model, X, y, cv=LeaveOneOut())
    loo_results[name] = evaluate(y, y_pred, name)

print_table(loo_results, "LOO-CV Results")

# ── 5. Leave-One-Series-Out CV ────────────────────────────────────────────────

print("\n--- Leave-One-Series-Out Cross-Validation ---")

series_map = {
    "25CN-NBOH (J Med Chem 2022)":                     "25CN-NBOH",
    "Heterocyclic Gq-biased (ACS Med Chem Lett 2023)":  "Heterocyclic",
    "literature (Pottie & Poulie 2023)":                "Pottie2023",
}
df["series"] = df["source"].map(series_map)
series_list  = df["series"].unique()

loso_results = {name: {"y_true": [], "y_pred": [], "series_held": []}
                for name in make_models()}

for held_out in series_list:
    train_mask = df["series"] != held_out
    test_mask  = df["series"] == held_out
    X_train, X_test = X[train_mask], X[test_mask]
    y_train, y_test = y[train_mask], y[test_mask]

    print(f"\n  Held out: {held_out} (n={test_mask.sum()}), "
          f"train: {train_mask.sum()} compounds")

    for name, model in make_models().items():
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        loso_results[name]["y_true"].extend(y_test.tolist())
        loso_results[name]["y_pred"].extend(y_pred.tolist())
        loso_results[name]["series_held"].extend([held_out] * len(y_test))
        mae_fold = mean_absolute_error(y_test, y_pred)
        print(f"    {name:<12} MAE={mae_fold:.3f}  "
              f"range predicted: [{y_pred.min():.2f}, {y_pred.max():.2f}]  "
              f"range actual: [{y_test.min():.2f}, {y_test.max():.2f}]")

# Aggregate LOSO metrics
loso_agg = {}
for name, res in loso_results.items():
    yt = np.array(res["y_true"])
    yp = np.array(res["y_pred"])
    loso_agg[name] = evaluate(yt, yp, name)
    loso_agg[name]["y_pred_all"] = yp
    loso_agg[name]["y_true_all"] = yt
    loso_agg[name]["series_held"] = res["series_held"]

print_table(loso_agg, "LOSO-CV Results (generalisation across chemical series)")

# ── 6. Compare LOO vs LOSO ────────────────────────────────────────────────────

print("\nLOO vs LOSO comparison (R²):")
print(f"{'Model':<12} {'LOO R²':>8} {'LOSO R²':>9} {'Drop':>8}")
print("-"*42)
for name in loo_results:
    loo_r2  = loo_results[name]["r2"]
    loso_r2 = loso_agg[name]["r2"]
    drop    = loo_r2 - loso_r2
    flag    = " *** likely confounded" if drop > 0.4 else ""
    print(f"{name:<12} {loo_r2:>8.3f} {loso_r2:>9.3f} {drop:>8.3f}{flag}")

# ── 7. Plots ──────────────────────────────────────────────────────────────────

source_colors = {
    "25CN-NBOH":   "#4C72B0",
    "Heterocyclic": "#DD8452",
    "Pottie2023":  "#55A868",
}

# LOO predicted vs actual
fig, axes = plt.subplots(2, 2, figsize=(11, 9))
for ax, (name, res) in zip(axes.flatten(), loo_results.items()):
    yp = res["y_pred"]
    for series, grp in df.groupby("series"):
        idx = grp.index
        ax.scatter(y[idx], yp[idx], alpha=0.8, s=40,
                   label=series, color=source_colors.get(series, "gray"))
    lim = max(abs(y.min()), abs(y.max()), abs(yp.min()), abs(yp.max())) * 1.1
    ax.plot([-lim, lim], [-lim, lim], "k--", lw=1, alpha=0.5)
    ax.axhline(0, color="gray", lw=0.5)
    ax.axvline(0, color="gray", lw=0.5)
    ax.set_xlabel("Actual log bias")
    ax.set_ylabel("Predicted log bias")
    ax.set_title(f"{name}  R²={res['r2']:.3f}  r={res['r']:.3f}")
    ax.legend(fontsize=6)
plt.suptitle("LOO-CV: predicted vs actual log bias", fontsize=11)
plt.tight_layout()
plt.savefig("analysis/qsar_loo_predictions.png", dpi=150)
print("\nSaved analysis/qsar_loo_predictions.png")

# LOSO predicted vs actual
fig2, axes2 = plt.subplots(2, 2, figsize=(11, 9))
for ax, (name, res) in zip(axes2.flatten(), loso_agg.items()):
    yt = res["y_true_all"]
    yp = res["y_pred_all"]
    sh = res["series_held"]
    for series in series_list:
        mask = np.array(sh) == series
        ax.scatter(np.array(yt)[mask], np.array(yp)[mask], alpha=0.8, s=40,
                   label=f"held: {series}", color=source_colors.get(series, "gray"))
    lim = max(abs(np.array(yt).min()), abs(np.array(yt).max()),
              abs(np.array(yp).min()), abs(np.array(yp).max())) * 1.1
    ax.plot([-lim, lim], [-lim, lim], "k--", lw=1, alpha=0.5)
    ax.axhline(0, color="gray", lw=0.5)
    ax.axvline(0, color="gray", lw=0.5)
    ax.set_xlabel("Actual log bias")
    ax.set_ylabel("Predicted log bias")
    ax.set_title(f"{name}  R²={res['r2']:.3f}  r={res['r']:.3f}")
    ax.legend(fontsize=6)
plt.suptitle("LOSO-CV: generalisation across chemical series", fontsize=11)
plt.tight_layout()
plt.savefig("analysis/qsar_loso_predictions.png", dpi=150)
print("Saved analysis/qsar_loso_predictions.png")

# Feature importance (RF, full dataset)
rf = RandomForestRegressor(n_estimators=200, max_features="sqrt",
                           random_state=42, n_jobs=-1)
rf.fit(X, y)
top_idx = np.argsort(rf.feature_importances_)[::-1][:20]
fig3, ax3 = plt.subplots(figsize=(10, 4))
ax3.bar(range(20), rf.feature_importances_[top_idx], color="#4C72B0", alpha=0.8)
ax3.set_xticks(range(20))
ax3.set_xticklabels([f"bit {i}" for i in top_idx], rotation=45, ha="right", fontsize=8)
ax3.set_ylabel("Feature importance")
ax3.set_title("Top 20 Morgan fingerprint bits (RF, full dataset)")
plt.tight_layout()
plt.savefig("analysis/feature_importance.png", dpi=150)
print("Saved analysis/feature_importance.png")

# ── 8. Report ─────────────────────────────────────────────────────────────────

best_loo  = max(loo_results,  key=lambda k: loo_results[k]["r2"])
best_loso = max(loso_agg,     key=lambda k: loso_agg[k]["r2"])

report = f"""
5-HT2A Biased Agonism QSAR — Results
======================================

Dataset: {len(df)} compounds, Morgan fingerprints (ECFP4 radius=2, {N_BITS} bits)
Informative bits (non-zero variance): {nonzero.sum()}
Target: log10 bias factor (continuous regression)

LOO-CV Results (within-dataset interpolation):
{'Model':<12} {'R²':>8} {'MAE':>8} {'Pearson r':>10} {'p-value':>14}
{'-'*56}
{''.join(f"{n:<12} {v['r2']:>8.3f} {v['mae']:>8.3f} {v['r']:>10.3f} {v['p']:>14.2e}{chr(10)}" for n, v in loo_results.items())}
Best LOO model: {best_loo}  R²={loo_results[best_loo]['r2']:.3f}

LOSO-CV Results (cross-series generalisation):
{'Model':<12} {'R²':>8} {'MAE':>8} {'Pearson r':>10} {'p-value':>14}
{'-'*56}
{''.join(f"{n:<12} {v['r2']:>8.3f} {v['mae']:>8.3f} {v['r']:>10.3f} {v['p']:>14.2e}{chr(10)}" for n, v in loso_agg.items())}
Best LOSO model: {best_loso}  R²={loso_agg[best_loso]['r2']:.3f}

LOO vs LOSO R² drop:
{'Model':<12} {'LOO R²':>8} {'LOSO R²':>9} {'Drop':>8}
{'-'*42}
{''.join(f"{n:<12} {loo_results[n]['r2']:>8.3f} {loso_agg[n]['r2']:>9.3f} {loo_results[n]['r2']-loso_agg[n]['r2']:>8.3f}{chr(10)}" for n in loo_results)}

Interpretation:
  A large LOO-LOSO drop (>0.4) suggests inter-series confounding:
  the model learns scaffold identity rather than the SAR within scaffolds.
  A small drop suggests the model captures genuine structure-bias signal
  that generalises across chemical series.

Chemical series in dataset:
{df.groupby('series')['log_bias'].describe()[['count','mean','std','min','max']].to_string()}

Notes:
  - p=0.00e+00 means p < 1e-300 (machine epsilon); not a calculation error
  - ChEMBL compounds use Emax-only bias; Pottie compounds use Emax/EC50
  - LOO-CV is appropriate for small datasets but optimistic vs true external test
"""

print(report)
with open("analysis/05_qsar_report.txt", "w") as f:
    f.write(report)
print("Done.")
