"""
Script 06: Within-series QSAR models with sensitivity analysis.

Trains and evaluates a separate QSAR model for each chemical series
using LOO-CV within that series. Also runs a sensitivity analysis
excluding known structural outliers from each series.

Structural outliers excluded in sensitivity analysis:
  25CN-NBOH series: SEROTONIN (endogenous neurotransmitter, structurally
    unlike the N-benzyl phenethylamine series; included in ChEMBL paper
    as a reference compound only)

Outputs:
  analysis/06_within_series_report.txt
  analysis/within_series_predictions.png
  analysis/within_series_sensitivity.png
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

# ── 1. Load and fingerprint ───────────────────────────────────────────────────

df = pd.read_csv("analysis/combined_dataset.csv")
df = df[df["log_bias"].notna() & df["smiles"].notna()].reset_index(drop=True)

series_map = {
    "25CN-NBOH (J Med Chem 2022)":                     "25CN-NBOH",
    "Heterocyclic Gq-biased (ACS Med Chem Lett 2023)":  "Heterocyclic",
    "literature (Pottie & Poulie 2023)":                "Pottie2023",
}
df["series"] = df["source"].map(series_map)

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

X_all = np.array(fps)
df    = df.loc[valid_idx].reset_index(drop=True)
y_all = df["log_bias"].values

print(f"Loaded {len(df)} compounds across {df['series'].nunique()} series")

# Structural outliers to exclude per series in sensitivity analysis
OUTLIERS = {
    "25CN-NBOH": {"SEROTONIN"},  # reference compound, not an NBOMe/NBOH
}

# ── 2. Models ─────────────────────────────────────────────────────────────────

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

# ── 3. Fit one series ─────────────────────────────────────────────────────────

def run_series(X_s, y_s, df_s, series, label=""):
    n = len(y_s)
    nonzero = X_s.var(axis=0) > 0
    X_nz    = X_s[:, nonzero]

    print(f"\n{'='*60}")
    print(f"Series: {series}{label}  (n={n}, {nonzero.sum()} informative bits)")
    print(f"Log bias range: {y_s.min():.3f} to {y_s.max():.3f}  std={y_s.std():.3f}")
    print(f"{'='*60}")

    model_results = {}
    for name, model in make_models().items():
        note = " (unreliable at n<15)" if n < 15 and name in ("RF", "GradBoost") else ""
        y_pred = cross_val_predict(model, X_nz, y_s, cv=LeaveOneOut())
        r2  = r2_score(y_s, y_pred)
        mae = mean_absolute_error(y_s, y_pred)
        r, p = stats.pearsonr(y_s, y_pred)
        model_results[name] = {
            "r2": r2, "mae": mae, "r": r, "p": p,
            "y_pred": y_pred, "note": note
        }
        print(f"  {name:<12} R²={r2:>7.3f}  MAE={mae:.3f}  "
              f"r={r:.3f}  p={p:.3e}{note}")

    best = max(("Ridge", "Lasso"),
               key=lambda k: model_results[k]["r2"])
    y_pred_best = model_results[best]["y_pred"]
    df_s = df_s.copy()
    df_s["y_pred"]   = y_pred_best
    df_s["residual"] = np.abs(y_s - y_pred_best)

    print(f"\n  Best linear model: {best}")
    print(f"  Worst predicted compounds:")
    worst = df_s.nlargest(3, "residual")[
        ["compound_name", "log_bias", "y_pred", "residual"]]
    print(worst.to_string(index=False))

    return {
        "models":  model_results,
        "df":      df_s,
        "y":       y_s,
        "best":    best,
        "n":       n,
        "nonzero": nonzero.sum(),
    }

# ── 4. Main run + sensitivity analysis ───────────────────────────────────────

series_results     = {}  # full series
sensitivity_results = {}  # outliers removed

for series in df["series"].unique():
    mask = df["series"] == series
    X_s  = X_all[mask]
    y_s  = y_all[mask]
    df_s = df[mask].reset_index(drop=True)

    # Main run
    series_results[series] = run_series(X_s, y_s, df_s, series)

    # Sensitivity: remove structural outliers if any defined for this series
    outliers = OUTLIERS.get(series, set())
    if outliers:
        keep = ~df_s["compound_name"].isin(outliers)
        n_removed = (~keep).sum()
        removed_names = df_s.loc[~keep, "compound_name"].tolist()
        print(f"\n  --- Sensitivity analysis: removing {n_removed} outlier(s): "
              f"{removed_names} ---")
        X_sens = X_s[keep]
        y_sens = y_s[keep]
        df_sens = df_s[keep].reset_index(drop=True)
        sensitivity_results[series] = run_series(
            X_sens, y_sens, df_sens, series,
            label=f" [excl. {', '.join(removed_names)}]")
        sensitivity_results[series]["removed"] = removed_names
    else:
        sensitivity_results[series] = None

# ── 5. Summary comparison table ───────────────────────────────────────────────

print("\n" + "="*70)
print("SUMMARY: Within-series Ridge R² (full vs sensitivity)")
print("="*70)
print(f"{'Series':<16} {'n_full':>6} {'R²_full':>9} "
      f"{'n_sens':>8} {'R²_sens':>9} {'Δ R²':>8}  Outliers removed")
print("-"*70)
for series in series_results:
    full = series_results[series]
    sens = sensitivity_results[series]
    r2_full = full["models"]["Ridge"]["r2"]
    if sens is not None:
        r2_sens = sens["models"]["Ridge"]["r2"]
        delta   = r2_sens - r2_full
        removed = ", ".join(sens["removed"])
        print(f"{series:<16} {full['n']:>6} {r2_full:>9.3f} "
              f"{sens['n']:>8} {r2_sens:>9.3f} {delta:>+8.3f}  {removed}")
    else:
        print(f"{series:<16} {full['n']:>6} {r2_full:>9.3f} "
              f"{'—':>8} {'—':>9} {'—':>8}  no outliers defined")

# ── 6. Plots ──────────────────────────────────────────────────────────────────

series_colors = {
    "25CN-NBOH":   "#4C72B0",
    "Heterocyclic": "#DD8452",
    "Pottie2023":  "#55A868",
}

def scatter_panel(ax, y_true, y_pred, df_s, series, title, color):
    ax.scatter(y_true, y_pred, color=color, alpha=0.85, s=55, zorder=3)
    for _, row in df_s.iterrows():
        ax.annotate(row["compound_name"][:10],
                    (row["log_bias"], row["y_pred"]),
                    fontsize=5, alpha=0.7,
                    xytext=(3, 3), textcoords="offset points")
    lim = max(abs(y_true.min()), abs(y_true.max()),
              abs(y_pred.min()), abs(y_pred.max())) * 1.15
    ax.plot([-lim, lim], [-lim, lim], "k--", lw=1, alpha=0.5)
    ax.axhline(0, color="gray", lw=0.5)
    ax.axvline(0, color="gray", lw=0.5)
    ax.set_xlabel("Actual log bias")
    ax.set_ylabel("Predicted log bias")
    ax.set_title(title)

# Main results
fig, axes = plt.subplots(1, 3, figsize=(14, 5))
for ax, (series, res) in zip(axes, series_results.items()):
    best   = res["best"]
    r2     = res["models"][best]["r2"]
    r      = res["models"][best]["r"]
    scatter_panel(ax, res["y"], res["models"][best]["y_pred"],
                  res["df"], series,
                  f"{series}\n{best} R²={r2:.3f} r={r:.3f} n={res['n']}",
                  series_colors.get(series, "gray"))
plt.suptitle("Within-series LOO-CV: predicted vs actual log bias", fontsize=11)
plt.tight_layout()
plt.savefig("analysis/within_series_predictions.png", dpi=150)
print("\nSaved analysis/within_series_predictions.png")

# Sensitivity analysis (only series with outliers removed)
sens_series = {s: r for s, r in sensitivity_results.items() if r is not None}
if sens_series:
    fig2, axes2 = plt.subplots(1, len(sens_series), figsize=(6*len(sens_series), 5))
    if len(sens_series) == 1:
        axes2 = [axes2]
    for ax, (series, res) in zip(axes2, sens_series.items()):
        best = res["best"]
        r2   = res["models"][best]["r2"]
        r    = res["models"][best]["r"]
        full_r2 = series_results[series]["models"][best]["r2"]
        scatter_panel(ax, res["y"], res["models"][best]["y_pred"],
                      res["df"], series,
                      f"{series} [excl. {', '.join(res['removed'])}]\n"
                      f"{best} R²={r2:.3f} (full: {full_r2:.3f}) n={res['n']}",
                      series_colors.get(series, "gray"))
    plt.suptitle("Sensitivity analysis: structural outliers removed", fontsize=11)
    plt.tight_layout()
    plt.savefig("analysis/within_series_sensitivity.png", dpi=150)
    print("Saved analysis/within_series_sensitivity.png")

# ── 7. Report ─────────────────────────────────────────────────────────────────

report_lines = ["""
Within-Series QSAR — LOO-CV Results with Sensitivity Analysis
==============================================================
Each model trained and evaluated within one chemical series only,
eliminating inter-series scaffold confounding.
"""]

for series in series_results:
    full = series_results[series]
    sens = sensitivity_results[series]
    report_lines.append(f"\nSeries: {series}  (n={full['n']}, "
                        f"{full['nonzero']} informative bits)")
    report_lines.append(f"  Log bias range: {full['y'].min():.3f} "
                        f"to {full['y'].max():.3f}")
    report_lines.append(f"\n  {'Model':<12} {'R²':>8} {'MAE':>8} "
                        f"{'Pearson r':>10} {'p-value':>14}")
    report_lines.append(f"  {'-'*56}")
    for name, r in full["models"].items():
        report_lines.append(
            f"  {name:<12} {r['r2']:>8.3f} {r['mae']:>8.3f} "
            f"{r['r']:>10.3f} {r['p']:>14.2e}{r['note']}")
    report_lines.append(f"\n  Best linear model: {full['best']}  "
                        f"R²={full['models'][full['best']]['r2']:.3f}")

    if sens is not None:
        report_lines.append(f"\n  Sensitivity analysis (excl. "
                            f"{', '.join(sens['removed'])}):")
        report_lines.append(f"  n={sens['n']}  "
                            f"{'Model':<12} {'R²':>8} {'MAE':>8} "
                            f"{'Pearson r':>10} {'p-value':>14}")
        report_lines.append(f"  {'-'*56}")
        for name, r in sens["models"].items():
            report_lines.append(
                f"  {' ':<14}{name:<12} {r['r2']:>8.3f} {r['mae']:>8.3f} "
                f"{r['r']:>10.3f} {r['p']:>14.2e}{r['note']}")
        delta = (sens["models"]["Ridge"]["r2"] -
                 full["models"]["Ridge"]["r2"])
        report_lines.append(f"\n  Ridge R² change after removing outlier: "
                            f"{delta:+.3f}")

report_lines.append("""
Interpretation:
  R² > 0.5 : meaningful within-series SAR signal
  R² 0.2-0.5: weak but present signal
  R² < 0.2 : insufficient signal at this n

Note: With n=12-17, LOO-CV R² estimates have high variance.
A ~0.1 difference between models is not meaningful.
The key question is whether R² is clearly above zero.
""")

report = "\n".join(report_lines)
print(report)
with open("analysis/06_within_series_report.txt", "w") as f:
    f.write(report)
print("Done.")
