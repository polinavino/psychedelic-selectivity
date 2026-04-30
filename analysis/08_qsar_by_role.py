import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score, KFold
import matplotlib.pyplot as plt

df = pd.read_csv('data/selectivity_with_roles.csv').dropna(subset=['smiles'])

def get_fingerprints(smiles_list):
    fps = []
    valid = []
    for i, smi in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smi)
        if mol:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024)
            fps.append(list(fp))
            valid.append(i)
    return np.array(fps), valid

cv = KFold(n_splits=5, shuffle=True, random_state=42)
rf = RandomForestRegressor(n_estimators=200, random_state=42, n_jobs=-1)

results = {}

for role in ['all', 'agonist', 'antagonist']:
    if role == 'all':
        subset = df
    else:
        subset = df[df['role'] == role]
    
    if len(subset) < 50:
        print(f"{role}: too few samples ({len(subset)})")
        continue
    
    X_fp, valid = get_fingerprints(subset['smiles'].tolist())
    subset = subset.iloc[valid].reset_index(drop=True)
    y = subset['selectivity_ratio'].values
    
    scores = cross_val_score(rf, X_fp, y, cv=cv, scoring='r2')
    results[role] = {'n': len(subset), 'r2': scores.mean(), 'std': scores.std()}
    print(f"{role:12s} (n={len(subset):3d}): R² = {scores.mean():.3f} ± {scores.std():.3f}")

# Plot selectivity distributions by role
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

colors = {'agonist': 'tomato', 'antagonist': 'steelblue'}
for role, color in colors.items():
    subset = df[df['role'] == role]['selectivity_ratio']
    axes[0].hist(subset, bins=30, alpha=0.6, color=color, label=f"{role} (n={len(subset)})")
axes[0].axvline(0, color='black', linestyle='--')
axes[0].set_xlabel('Selectivity (pKi_2A - pKi_2B)')
axes[0].set_ylabel('Count')
axes[0].set_title('Selectivity distribution by pharmacological role')
axes[0].legend()

# R² comparison
roles = list(results.keys())
r2s = [results[r]['r2'] for r in roles]
stds = [results[r]['std'] for r in roles]
axes[1].bar(roles, r2s, yerr=stds, color=['grey','tomato','steelblue'], alpha=0.8, capsize=5)
axes[1].set_ylabel('Cross-validated R²')
axes[1].set_title('QSAR model performance by role')
axes[1].axhline(0, color='black', linewidth=0.5)

plt.tight_layout()
plt.savefig('analysis/qsar_by_role.png', dpi=150, bbox_inches='tight')
plt.close()
print("\nSaved analysis/qsar_by_role.png")
