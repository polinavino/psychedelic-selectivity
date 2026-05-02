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

## Serotonin receptor selectivity framework (Option A)

**Script:** `analysis/09_full_serotonin_family.py`, `analysis/10_selectivity_framework.py`  
**Data:** ChEMBL Ki values for 13 human serotonin receptors (5HT3B excluded — no data)

This analysis applies the selectivity definition framework from the companion
kinase selectivity paper to the full serotonin receptor family, asking whether
the same definitional instabilities appear in receptor pharmacology.

### Data

Ki data downloaded from ChEMBL for all human serotonin receptors:

| Receptor | ChEMBL ID | Compounds |
|----------|-----------|-----------|
| 5-HT1A | CHEMBL214 | 4,859 |
| 5-HT1B | CHEMBL1898 | 751 |
| 5-HT1D | CHEMBL1983 | 916 |
| 5-HT1E | CHEMBL2182 | 72 |
| 5-HT1F | CHEMBL1805 | 113 |
| 5-HT2A | CHEMBL224 | 4,602 |
| 5-HT2B | CHEMBL1833 | 1,554 |
| 5-HT2C | CHEMBL225 | 2,526 |
| 5-HT3A | CHEMBL1899 | 579 |
| 5-HT4 | CHEMBL1875 | 424 |
| 5-HT5A | CHEMBL3426 | 344 |
| 5-HT6 | CHEMBL3371 | 3,782 |
| 5-HT7 | CHEMBL3155 | 2,633 |

Total: 13,584 unique compounds. 297 compounds have Ki values at 5-HT2A,
5-HT2B, and at least 2 other receptors — used for selectivity analysis.

### Selectivity scores for psilocin

| Metric | Value | Interpretation |
|--------|-------|----------------|
| S-score | 1.000 | All 10 tested receptors above threshold |
| Entropy | 3.288 bits | High — distributed response |
| Gini | 0.118 | Low — binding spread uniformly |
| Ratio | 1.080 | Top receptor barely stronger than second |

Psilocin binds all tested serotonin receptors with similar affinity (pKi
6.5-8.3). Its highest affinity is for 5-HT2B (pKi=8.34) and 5-HT1D
(pKi=7.72), not 5-HT2A (pKi=6.72). By all four definitions, psilocin is
a non-selective serotonergic compound — inconsistent with its common
description as a "5-HT2A agonist."

### Correlations between selectivity definitions

| | S-score | Entropy | Gini | Ratio |
|-|---------|---------|------|-------|
| S-score | 1.000 | -0.063 | -0.682 | -0.158 |
| Entropy | -0.063 | 1.000 | 0.100 | -0.047 |
| Gini | -0.682 | 0.100 | 1.000 | 0.313 |
| Ratio | -0.158 | -0.047 | 0.313 | 1.000 |

### Key finding

Entropy is uncorrelated with all other definitions (r ≈ 0), while S-score
and Gini are strongly anti-correlated (r = -0.682). This reproduces the
two-cluster structure identified in the kinase selectivity paper — one cluster
of mutually consistent definitions (S-score, Gini, ratio) and entropy as an
outlier measuring something categorically different.

This finding generalizes the kinase paper result to serotonin receptor
pharmacology, suggesting that entropy's divergence from other selectivity
definitions is a property of the mathematical definition, not an artifact
of kinase binding data specifically.
