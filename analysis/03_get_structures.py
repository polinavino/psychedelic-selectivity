import pandas as pd
from chembl_webresource_client.new_client import new_client
import time

both = pd.read_csv('data/selectivity_data.csv')
molecule = new_client.molecule

print(f"Getting SMILES for {len(both)} compounds...")

smiles_list = []
batch_size = 50

for i in range(0, len(both), batch_size):
    batch = both['molecule_chembl_id'].iloc[i:i+batch_size].tolist()
    try:
        mols = molecule.filter(molecule_chembl_id__in=batch).only([
            'molecule_chembl_id',
            'pref_name',
            'molecule_structures'
        ])
        for m in mols:
            smiles = None
            if m.get('molecule_structures'):
                smiles = m['molecule_structures'].get('canonical_smiles')
            smiles_list.append({
                'molecule_chembl_id': m['molecule_chembl_id'],
                'pref_name': m.get('pref_name'),
                'smiles': smiles
            })
    except Exception as e:
        print(f"Error at batch {i}: {e}")
    
    if (i // batch_size) % 5 == 0:
        print(f"  {i+batch_size}/{len(both)} done...")
    time.sleep(0.2)

smiles_df = pd.DataFrame(smiles_list)
result = both.merge(smiles_df[['molecule_chembl_id','smiles']], 
                    on='molecule_chembl_id', how='left')

print(f"Compounds with SMILES: {result['smiles'].notna().sum()}")
result.to_csv('data/selectivity_with_smiles.csv', index=False)
print("Saved data/selectivity_with_smiles.csv")
