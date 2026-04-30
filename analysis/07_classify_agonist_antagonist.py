import pandas as pd
import numpy as np

ht2a_func = pd.read_csv('data/ht2a_functional_raw.csv')
ht2b_func = pd.read_csv('data/ht2b_functional_raw.csv')

def classify_from_description(desc):
    if pd.isna(desc):
        return None
    desc = desc.lower()
    if 'antagonist' in desc:
        return 'antagonist'
    elif 'agonist' in desc:
        return 'agonist'
    return None

ht2a_func['role'] = ht2a_func['assay_description'].apply(classify_from_description)
ht2b_func['role'] = ht2b_func['assay_description'].apply(classify_from_description)

print("5-HT2A functional role distribution:")
print(ht2a_func['role'].value_counts(dropna=False))

print("\n5-HT2B functional role distribution:")
print(ht2b_func['role'].value_counts(dropna=False))

# Get compounds with clear agonist/antagonist labels
ht2a_labeled = ht2a_func.dropna(subset=['role', 'molecule_chembl_id'])
ht2b_labeled = ht2b_func.dropna(subset=['role', 'molecule_chembl_id'])

# Take majority vote per compound
ht2a_role = ht2a_labeled.groupby('molecule_chembl_id')['role'].agg(
    lambda x: x.value_counts().index[0]).reset_index()
ht2a_role.columns = ['molecule_chembl_id', 'role_2a']

ht2b_role = ht2b_labeled.groupby('molecule_chembl_id')['role'].agg(
    lambda x: x.value_counts().index[0]).reset_index()
ht2b_role.columns = ['molecule_chembl_id', 'role_2b']

print(f"\nUnique compounds with 5-HT2A role: {len(ht2a_role)}")
print(f"Unique compounds with 5-HT2B role: {len(ht2b_role)}")

# Merge with selectivity data
both = pd.read_csv('data/selectivity_with_descriptors.csv')
both = both.merge(ht2a_role, on='molecule_chembl_id', how='left')
both = both.merge(ht2b_role, on='molecule_chembl_id', how='left')

# Consistent role across both receptors
both['role'] = np.where(
    both['role_2a'] == both['role_2b'], both['role_2a'],
    np.where(both['role_2a'].notna(), both['role_2a'],
    np.where(both['role_2b'].notna(), both['role_2b'], None))
)

print("\nRole distribution in selectivity dataset:")
print(both['role'].value_counts(dropna=False))

both.to_csv('data/selectivity_with_roles.csv', index=False)
print("\nSaved data/selectivity_with_roles.csv")

# Compare selectivity by role
print("\nSelectivity ratio by role:")
print(both.groupby('role')['selectivity_ratio'].describe().round(3))
