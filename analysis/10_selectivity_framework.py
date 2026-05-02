import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

matrix = pd.read_csv('data/serotonin_receptor_matrix.csv')

receptor_cols = [c for c in matrix.columns if c != 'molecule_chembl_id']
print(f"Receptors: {receptor_cols}")
print(f"Total compounds: {len(matrix)}")

# ── Focus on compounds with 2A + 2B + at least 3 others ─────────────────────
has_2a_2b = matrix[['5HT2A','5HT2B']].notna().all(axis=1)
n_receptors = matrix[receptor_cols].notna().sum(axis=1)
subset = matrix[has_2a_2b & (n_receptors >= 4)].copy()
print(f"Compounds with 2A+2B+>=2 others: {len(subset)}")

# ── Selectivity definitions (analog of kinase paper) ─────────────────────────
# Focus on 5HT2A as primary target, all others as off-targets
# Higher pKi at 2A, lower at others = more selective

def s_score(row, primary='5HT2A', threshold=6.0):
    """Fraction of tested receptors above threshold (lower = more selective)"""
    vals = row[receptor_cols].dropna()
    if len(vals) == 0:
        return np.nan
    return (vals > threshold).sum() / len(vals)

def entropy_score(row, primary='5HT2A', beta=5.0):
    """Shannon entropy of normalized pKi distribution (lower = more selective)"""
    vals = np.array(row[receptor_cols].dropna().values, dtype=float)
    shifted = np.maximum(vals - beta, 0)
    total = shifted.sum()
    if total == 0:
        return 0.0
    p = shifted / total
    p = p[p > 0]
    return float(-np.sum(p * np.log2(p)))

def gini_score(row):
    """Gini coefficient (higher = more selective)"""
    vals = row[receptor_cols].dropna().values
    vals = np.maximum(vals - 5.0, 0)
    if vals.sum() == 0:
        return 0.0
    n = len(vals)
    sorted_vals = np.sort(vals)
    indices = np.arange(1, n + 1)
    return (2 * np.sum(indices * sorted_vals)) / (n * np.sum(sorted_vals)) - (n + 1) / n

def ratio_score(row, primary='5HT2A'):
    """Ratio of top target to second target pKi"""
    vals = row[receptor_cols].dropna().sort_values(ascending=False)
    if len(vals) < 2 or vals.iloc[1] <= 5.0:
        return np.nan
    return vals.iloc[0] / vals.iloc[1]

subset['s_score']       = subset.apply(s_score, axis=1)
subset['entropy']       = subset.apply(entropy_score, axis=1)
subset['gini']          = subset.apply(gini_score, axis=1)
subset['ratio']         = subset.apply(ratio_score, axis=1)
subset['n_receptors']   = subset[receptor_cols].notna().sum(axis=1)

# ── Check known psychedelics ──────────────────────────────────────────────────
known = {
    'LSD':        'CHEMBL463207',
    'Psilocin':   'CHEMBL65547',
    'Mescaline':  'CHEMBL26687',
    'MDMA':       'CHEMBL43048',
    'Ketanserin': 'CHEMBL279218',
    'Risperidone':'CHEMBL1721',
}

print("\n── Known compounds in matrix ──")
for name, cid in known.items():
    row = subset[subset['molecule_chembl_id'] == cid]
    if len(row) > 0:
        r = row.iloc[0]
        receptors_tested = [c for c in receptor_cols if not pd.isna(r[c])]
        print(f"\n{name}:")
        print(f"  Receptors tested: {receptors_tested}")
        for rc in receptors_tested:
            print(f"    {rc}: {r[rc]:.2f}")
        print(f"  S-score:  {r['s_score']:.3f}")
        print(f"  Entropy:  {r['entropy']:.3f}")
        print(f"  Gini:     {r['gini']:.3f}")
        print(f"  Ratio:    {r['ratio']:.3f}")
    else:
        print(f"{name}: not found in subset")

# ── Correlation between definitions ──────────────────────────────────────────
print("\n── Correlations between selectivity definitions ──")
score_cols = ['s_score', 'entropy', 'gini', 'ratio']
valid = subset[score_cols].dropna()
print(f"Compounds with all 4 scores: {len(valid)}")
corr = valid.corr(method='spearman')
print(corr.round(3))

subset.to_csv('data/selectivity_framework_results.csv', index=False)
print("\nSaved data/selectivity_framework_results.csv")

# ── Plot correlation matrix ───────────────────────────────────────────────────
import matplotlib.pyplot as plt
import seaborn as sns

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Correlation heatmap
sns.heatmap(corr, annot=True, fmt='.3f', cmap='RdBu_r', center=0,
            ax=axes[0], square=True, vmin=-1, vmax=1)
axes[0].set_title('Spearman correlations between\nselectivity definitions (serotonin receptors)')

# Scatter: entropy vs S-score colored by Gini
sc = axes[1].scatter(valid['s_score'], valid['entropy'],
                     c=valid['gini'], cmap='RdYlBu', alpha=0.5, s=20)
plt.colorbar(sc, ax=axes[1], label='Gini coefficient')
axes[1].set_xlabel('S-score')
axes[1].set_ylabel('Entropy (bits)')
axes[1].set_title('Entropy vs S-score\n(colored by Gini)')

plt.tight_layout()
plt.savefig('analysis/selectivity_framework.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved analysis/selectivity_framework.png")
