from chembl_webresource_client.new_client import new_client
import pandas as pd

activity = new_client.activity

def get_functional_data(target_id, target_name):
    print(f"Fetching functional data for {target_name}...")
    data = activity.filter(
        target_chembl_id=target_id,
        assay_type='F',  # Functional assays
    ).only([
        'molecule_chembl_id',
        'molecule_pref_name',
        'standard_type',
        'standard_value',
        'standard_units',
        'pchembl_value',
        'assay_description',
        'assay_chembl_id',
    ])
    df = pd.DataFrame(data)
    df['target'] = target_name
    print(f"  {len(df)} records")
    print(f"  Standard types: {df['standard_type'].value_counts().head(10).to_dict()}")
    return df

ht2a_func = get_functional_data('CHEMBL224', '5HT2A')
ht2b_func = get_functional_data('CHEMBL1833', '5HT2B')

ht2a_func.to_csv('data/ht2a_functional_raw.csv', index=False)
ht2b_func.to_csv('data/ht2b_functional_raw.csv', index=False)

print("\nSample assay descriptions for 5-HT2A:")
print(ht2a_func['assay_description'].dropna().head(10).tolist())
