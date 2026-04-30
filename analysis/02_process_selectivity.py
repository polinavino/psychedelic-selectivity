import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import matplotlib.pyplot as plt

# ── 1. Load and clean ─────────────────────────────────────────────────────────
ht2a = pd.read_csv('data/ht2a_ki_raw.csv')
ht2b = pd.read_csv('data/ht2b_ki_raw.csv')

# Keep only valid pChEMBL values (pKi = -log10(Ki in M))
ht2a = ht2a.dropna(subset=['pchembl_value', 'molecule_chembl_id'])
ht2b = ht2b.dropna(subset=['pchembl_value', 'molecule_chembl_id'])

ht2a['pchembl_value'] = pd.to_numeric(ht2a['pchembl_value'], errors='coerce')
ht2b['pchembl_value'] = pd.to_numeric(ht2b['pchembl_value'], errors='coerce')

ht2a = ht2a.dropna(subset=['pchembl_value'])
ht2b = ht2b.dropna(subset=['pchembl_value'])

# Take median pKi per compound (multiple assays)
ht2a_med = ht2a.groupby('molecule_chembl_id').agg(
    pki_2a = ('pchembl_value', 'median'),
    n_assays_2a = ('pchembl_value', 'count'),
    name = ('molecule_pref_name', 'first')
).reset_index()

ht2b_med = ht2b.groupby('molecule_chembl_id').agg(
    pki_2b = ('pchembl_value', 'median'),
    n_assays_2b = ('pchembl_value', 'count')
).reset_index()

print(f"Unique compounds with 5-HT2A Ki: {len(ht2a_med)}")
print(f"Unique compounds with 5-HT2B Ki: {len(ht2b_med)}")

# ── 2. Find compounds with both measurements ───────────────────────────────────
both = ht2a_med.merge(ht2b_med, on='molecule_chembl_id', how='inner')
print(f"Compounds with both 5-HT2A and 5-HT2B Ki: {len(both)}")

# ── 3. Compute selectivity ────────────────────────────────────────────────────
# Selectivity for 5-HT2A over 5-HT2B
# Higher = more selective for 5-HT2A (safer, less cardiac risk)
# pKi_2A - pKi_2B > 0 means more potent at 2A than 2B

both['selectivity_ratio'] = both['pki_2a'] - both['pki_2b']  # log ratio

print(f"\nSelectivity ratio (pKi_2A - pKi_2B):")
print(f"  Mean: {both['selectivity_ratio'].mean():.3f}")
print(f"  Std:  {both['selectivity_ratio'].std():.3f}")
print(f"  Min:  {both['selectivity_ratio'].min():.3f}")
print(f"  Max:  {both['selectivity_ratio'].max():.3f}")

print(f"\nCompounds more potent at 2A than 2B (selective): {(both['selectivity_ratio'] > 0).sum()}")
print(f"Compounds more potent at 2B than 2A (risky):    {(both['selectivity_ratio'] < 0).sum()}")

# ── 4. Check known psychedelics ───────────────────────────────────────────────
known = {
    'LSD': 'CHEMBL600',
    'Psilocin': 'CHEMBL6895',
    'DMT': 'CHEMBL14113',
    'Mescaline': 'CHEMBL274978',
    'DOI': 'CHEMBL260748',
    'MDMA': 'CHEMBL400',
}

print("\n── Known psychedelics ──")
for name, chembl_id in known.items():
    row = both[both['molecule_chembl_id'] == chembl_id]
    if len(row) > 0:
        r = row.iloc[0]
        print(f"{name:12s}: pKi_2A={r['pki_2a']:.2f}, pKi_2B={r['pki_2b']:.2f}, selectivity={r['selectivity_ratio']:+.2f}")
    else:
        print(f"{name:12s}: not found in both datasets")

# ── 5. Save ───────────────────────────────────────────────────────────────────
both.to_csv('data/selectivity_data.csv', index=False)
print(f"\nSaved data/selectivity_data.csv")

# ── 6. Plot ───────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

axes[0].scatter(both['pki_2a'], both['pki_2b'], alpha=0.3, s=10, color='steelblue')
axes[0].plot([4, 11], [4, 11], 'r--', alpha=0.5, label='No selectivity')
axes[0].set_xlabel('pKi 5-HT2A')
axes[0].set_ylabel('pKi 5-HT2B')
r, p = spearmanr(both['pki_2a'], both['pki_2b'])
axes[0].set_title(f'5-HT2A vs 5-HT2B binding\nSpearman r={r:.3f}, p={p:.2e}')
axes[0].legend()

axes[1].hist(both['selectivity_ratio'], bins=50, color='steelblue', alpha=0.7)
axes[1].axvline(0, color='red', linestyle='--', label='No selectivity')
axes[1].axvline(1, color='orange', linestyle=':', label='10x selective')
axes[1].axvline(-1, color='orange', linestyle=':')
axes[1].set_xlabel('Selectivity (pKi_2A - pKi_2B)')
axes[1].set_ylabel('Number of compounds')
axes[1].set_title('Distribution of 5-HT2A/2B selectivity')
axes[1].legend()

plt.tight_layout()
plt.savefig('analysis/selectivity_overview.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved analysis/selectivity_overview.png")
