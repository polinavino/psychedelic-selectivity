# Psychedelic Receptor Selectivity Analysis

QSAR analysis of 5-HT2A vs 5-HT2B serotonin receptor selectivity using
ChEMBL binding data, motivated by cardiac safety concerns in psychedelic
drug discovery.

## Background

Serotonergic psychedelics act primarily through the 5-HT2A receptor.
However, many compounds in this chemical space also bind the 5-HT2B receptor,
whose sustained activation causes cardiac valvulopathy (heart valve damage),
as seen with fenfluramine. Designing compounds with high 5-HT2A/5-HT2B
selectivity is a key safety objective in psychedelic drug discovery.

## Data

Binding affinity (Ki) data downloaded from ChEMBL:
- Human 5-HT2A (CHEMBL224): 6,170 measurements, 4,602 unique compounds
- Human 5-HT2B (CHEMBL1833): 1,785 measurements, 1,554 unique compounds
- Compounds with both measurements: 917

Functional assay data used to classify compounds as agonists or antagonists
based on assay description text parsing (majority vote per compound).

## Key findings

**Selectivity distribution:**
- Mean selectivity (pKi_2A - pKi_2B) = -0.26 across 917 compounds
- 576/917 compounds are more potent at 5-HT2B than 5-HT2A
- Only a minority achieve >10x selectivity for 5-HT2A

**Role-stratified selectivity:**
- Agonists (n=162): mean selectivity = -0.341 (favor 5-HT2B)
- Antagonists (n=87): mean selectivity = +0.089 (slight 5-HT2A preference)
- Agonists preferentially bind the cardiotoxic receptor — relevant safety signal

**QSAR model performance (5-fold CV R²):**

| Dataset | Features | R² |
|---------|----------|-----|
| All compounds | Physicochemical descriptors | 0.225 |
| All compounds | Morgan fingerprints | 0.462 |
| Agonists only | Morgan fingerprints | 0.180 |
| Antagonists only | Morgan fingerprints | 0.498 |

Structural features explain antagonist selectivity moderately well but
agonist selectivity poorly. Agonist 5-HT2A/2B selectivity likely depends
on functional selectivity (biased agonism) and receptor conformation not
captured by 2D fingerprints.

## Known psychedelics

| Compound | pKi_2A | pKi_2B | Selectivity |
|----------|--------|--------|-------------|
| LSD (Lysergide) | 8.50 | 7.52 | +0.98 |
| Psilocin | 6.72 | 8.34 | -1.62 |

LSD shows modest 5-HT2A preference. Psilocin's strong 5-HT2B preference
(~40x more potent at 2B than 2A) is notable given its therapeutic use, though
this may reflect assay variability across the limited measurements available.

## Limitations

- Dataset mixes agonists, antagonists, and partial agonists from heterogeneous assays
- Only 249/917 compounds have role labels from functional data
- 2D fingerprints cannot capture biased agonism or receptor conformation effects
- Classic psychedelics (DMT, mescaline, DOB series) have sparse ChEMBL coverage —
  PDSP Ki database would provide better coverage but requires manual download

## Next steps

- Add PDSP Ki database data for better psychedelic coverage
- Incorporate Emax values as a continuous agonism measure rather than binary classification
- Use 3D descriptors or docking scores from published 5-HT2A/2B crystal structures
- Apply the selectivity definition framework from the companion kinase selectivity
  paper to the receptor binding profiles across the full serotonin receptor family

## Scripts

1. `analysis/01_download_chembl.py` — fetch Ki data from ChEMBL API
2. `analysis/02_process_selectivity.py` — compute selectivity ratios
3. `analysis/03_get_structures.py` — fetch SMILES from ChEMBL
4. `analysis/04_qsar_model.py` — QSAR with physicochemical descriptors
5. `analysis/05_fingerprint_qsar.py` — QSAR with Morgan fingerprints
6. `analysis/06_get_functional_data.py` — fetch functional assay data
7. `analysis/07_classify_agonist_antagonist.py` — classify by role
8. `analysis/08_qsar_by_role.py` — role-stratified QSAR

## Connection to companion work

This project applies the receptor selectivity analysis framework to serotonin
receptors, complementing the kinase inhibitor selectivity paper
(ChemRxiv: https://doi.org/10.26434/chemrxiv.15001618/v1) which developed
formal desiderata for binding-based selectivity metrics. A natural extension
would apply those desiderata to the full serotonin receptor family binding
profiles of psychedelic compounds.
