from chembl_webresource_client.new_client import new_client
import pandas as pd

activity = new_client.activity

# Human 5-HT2A — CHEMBL224
# Human 5-HT2B — CHEMBL1833

def get_ki_data(target_id, target_name):
    print(f"Fetching Ki data for {target_name} ({target_id})...")
    data = activity.filter(
        target_chembl_id=target_id,
        standard_type='Ki',
        relation='=',
        assay_type='B'
    ).only([
        'molecule_chembl_id',
        'molecule_pref_name',
        'standard_value',
        'standard_units',
        'standard_type',
        'pchembl_value',
        'assay_chembl_id',
        'document_chembl_id',
    ])
    df = pd.DataFrame(data)
    df['target'] = target_name
    df['target_id'] = target_id
    print(f"  {len(df)} records")
    return df

ht2a = get_ki_data('CHEMBL224', '5HT2A')
ht2b = get_ki_data('CHEMBL1833', '5HT2B')

ht2a.to_csv('data/ht2a_ki_raw.csv', index=False)
ht2b.to_csv('data/ht2b_ki_raw.csv', index=False)

print(f"\nHT2A: {len(ht2a)} Ki measurements")
print(f"HT2B: {len(ht2b)} Ki measurements")
print(f"\nHT2A columns: {ht2a.columns.tolist()}")
print(f"\nHT2A sample:\n{ht2a[['molecule_pref_name','standard_value','pchembl_value']].head(10)}")
