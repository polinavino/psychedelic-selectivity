import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score, cross_val_predict, KFold
import matplotlib.pyplot as plt

df = pd.read_csv('data/selectivity_with_descriptors.csv').dropna(subset=['smiles'])

# Morgan fingerprints (radius=2, 1024 bits)
fps = []
valid_idx = []
for i, row in df.iterrows():
    mol = Chem.MolFromSmiles(row['smiles'])
    if mol:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024)
        fps.append(list(fp))
        valid_idx.append(i)

df = df.loc[valid_idx].reset_index(drop=True)
X_fp = np.array(fps)
y = df['selectivity_ratio'].values

print(f"Compounds with valid fingerprints: {len(df)}")

cv = KFold(n_splits=5, shuffle=True, random_state=42)

rf = RandomForestRegressor(n_estimators=200, random_state=42, n_jobs=-1)
scores = cross_val_score(rf, X_fp, y, cv=cv, scoring='r2')
print(f"Random Forest (Morgan FP): R² = {scores.mean():.3f} ± {scores.std():.3f}")

y_pred = cross_val_predict(rf, X_fp, y, cv=cv)

# Combined descriptors + fingerprints
desc_cols = ['MW', 'LogP', 'HBD', 'HBA', 'TPSA', 'RotBonds', 'ArRings', 'RingCount']
X_combined = np.hstack([X_fp, df[desc_cols].values])
scores_comb = cross_val_score(rf, X_combined, y, cv=cv, scoring='r2')
print(f"Random Forest (FP + desc): R² = {scores_comb.mean():.3f} ± {scores_comb.std():.3f}")

fig, ax = plt.subplots(figsize=(6, 5))
ax.scatter(y, y_pred, alpha=0.3, s=10, color='steelblue')
ax.plot([-3, 3], [-3, 3], 'r--', alpha=0.5)
ax.set_xlabel('Actual selectivity (pKi_2A - pKi_2B)')
ax.set_ylabel('Predicted selectivity')
ax.set_title(f'Random Forest (Morgan FP): R² = {scores.mean():.3f}')
plt.tight_layout()
plt.savefig('analysis/qsar_fingerprint.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved analysis/qsar_fingerprint.png")
