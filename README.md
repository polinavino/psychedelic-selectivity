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

* Human 5-HT2A (CHEMBL224): 6,170 measurements, 4,602 unique compounds
* Human 5-HT2B (CHEMBL1833): 1,785 measurements, 1,554 unique compounds
* Compounds with both measurements: 917

Functional assay data used to classify compounds as agonists or antagonists
based on assay description text parsing (majority vote per compound).

## Key findings

**Selectivity distribution:**

* Mean selectivity (pKi\_2A - pKi\_2B) = -0.26 across 917 compounds
* 576/917 compounds are more potent at 5-HT2B than 5-HT2A
* Only a minority achieve >10x selectivity for 5-HT2A

**Role-stratified selectivity:**

* Agonists (n=162): mean selectivity = -0.341 (favor 5-HT2B)
* Antagonists (n=87): mean selectivity = +0.089 (slight 5-HT2A preference)
* Agonists preferentially bind the cardiotoxic receptor — relevant safety signal

**QSAR model performance (5-fold CV R²):**

| Dataset | Features | R² |
| --- | --- | --- |
| All compounds | Physicochemical descriptors | 0.225 |
| All compounds | Morgan fingerprints | 0.462 |
| Agonists only | Morgan fingerprints | 0.180 |
| Antagonists only | Morgan fingerprints | 0.498 |

Structural features explain antagonist selectivity moderately well but
agonist selectivity poorly. Agonist 5-HT2A/2B selectivity likely depends
on functional selectivity (biased agonism) and receptor conformation not
captured by 2D fingerprints.

## Known psychedelics

| Compound | pKi\_2A | pKi\_2B | Selectivity |
| --- | --- | --- | --- |
| LSD (Lysergide) | 8.50 | 7.52 | +0.98 |
| Psilocin | 6.72 | 8.34 | -1.62 |

LSD shows modest 5-HT2A preference. Psilocin's strong 5-HT2B preference
(~40x more potent at 2B than 2A) is notable given its therapeutic use, though
this may reflect assay variability across the limited measurements available.

## Limitations

* Dataset mixes agonists, antagonists, and partial agonists from heterogeneous assays
* Only 249/917 compounds have role labels from functional data
* 2D fingerprints cannot capture biased agonism or receptor conformation effects
* Classic psychedelics (DMT, mescaline, DOB series) have sparse ChEMBL coverage —
  PDSP Ki database would provide better coverage but requires manual download

## Next steps

* Add PDSP Ki database data for better psychedelic coverage
* Incorporate Emax values as a continuous agonism measure rather than binary classification
* Use 3D descriptors or docking scores from published 5-HT2A/2B crystal structures
* Apply the selectivity definition framework from the companion kinase selectivity
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
(ChemRxiv: <https://doi.org/10.26434/chemrxiv.15001618/v1>) which developed
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
| --- | --- | --- |
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
| --- | --- | --- |
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

|  | S-score | Entropy | Gini | Ratio |
| --- | --- | --- | --- | --- |
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

---

## Biased agonism QSAR at 5-HT2A

> **In plain terms:** At the 5-HT2A receptor, the same compound can trigger
> two completely different cellular responses depending on which internal
> signalling pathway it activates. One pathway (Gq protein) is associated
> with hallucinations; the other (β-arrestin) is associated with
> neuroplasticity and antidepressant effects. Compounds that preferentially
> activate one pathway over the other are called *biased agonists*. This
> analysis asks whether a molecule's 2D structure is enough to predict which
> pathway it will prefer — a question with direct implications for designing
> non-hallucinogenic psychedelic-inspired therapeutics.

### Scientific motivation

Non-hallucinogenic psychedelics are an active area of drug discovery. At
the 5-HT2A receptor, agonists can activate two distinct intracellular
signalling pathways:

- **Gq protein signalling** — associated with hallucinogenic effects
- **β-arrestin-2 recruitment** — associated with neuroplasticity and
  potentially therapeutic effects independent of hallucination

Compounds that preferentially activate one pathway are called *biased
agonists*. Predicting bias from molecular structure would directly support
the rational design of compounds with the therapeutic profile of
psychedelics but without the hallucinogenic effects.

### Bias factor definition

For each compound tested in both a Gq-pathway assay and a β-arrestin assay,
we compute:

```
log_bias = log10( (Emax_βarr / EC50_βarr) / (Emax_Gq / EC50_Gq) )
```

- **log_bias > 0**: β-arrestin biased (potentially non-hallucinogenic)
- **log_bias < 0**: Gq biased (potentially hallucinogenic)
- **log_bias ≈ 0**: balanced (activates both pathways equally)

This is the transduction ratio formulation from the operational model of
agonism (Kenakin et al.). When EC50 is unavailable for a compound, we fall
back to the Emax-only ratio `log10(Emax_βarr / Emax_Gq)`.

### Data curation

A key methodological decision is that bias factors are only computed for
compounds measured within the same paper, using the same assay platform and
the same reference agonist. Cross-lab Emax comparisons are not made because
assay platforms (calcium flux, IP1 accumulation, BRET) are not
directly comparable even when reporting nominally identical quantities.

**ChEMBL core dataset** (scripts 01–03): All functional assay records for
CHEMBL224 (5-HT2A) were fetched from the ChEMBL API and classified by
pathway using keyword matching on assay descriptions. Of 295 unique functional
assays, 7 measured β-arrestin recruitment and 200 measured Gq-pathway
activation. Only 2 papers contained both pathway types, yielding 29 compounds
with paired Emax measurements — all using BRET assays in HEK293T cells, the
gold-standard platform for biased agonism measurement.

**Literature curation** (script 04): Table 1 from Pottie & Poulie et al.
(ACS Chem. Neurosci. 2023, DOI: 10.1021/acschemneuro.3c00267) was manually
transcribed, providing 17 additional compounds (phenethylamines and
N-benzyl derivatives, LSD as reference, all using NanoBiT miniGαq and
β-arrestin-2 recruitment assays).

**Combined dataset**: 46 compounds with SMILES, spanning three chemical
series.

| Source | n | Bias metric | Log bias range |
|---|---|---|---|
| 25CN-NBOH series (J Med Chem 2022) | 12 | Emax ratio | −0.30 to +0.67 |
| Heterocyclic Gq-biased (ACS Med Chem Lett 2023) | 17 | Emax ratio | −0.60 to +0.04 |
| Pottie & Poulie 2023 (phenethylamines) | 17 | Emax/EC50 ratio | +0.19 to +0.86 |

Note that the Pottie compounds use the full transduction ratio while the
ChEMBL compounds use Emax-only, because the two ChEMBL source papers used
different reference agonists (LSD vs serotonin), making EC50-based
transduction ratios incomparable across papers without re-normalisation.
This inconsistency in bias metric is a limitation of the dataset and adds
noise to any QSAR model trained on the combined data.

### QSAR results

Morgan fingerprints (ECFP4, radius=2, 2048 bits) were computed for all
46 compounds using RDKit. After removing zero-variance bits, 306 informative
bits remained. Four models were evaluated: Ridge regression, Lasso, Random
Forest, and Gradient Boosting.

**Cross-validation scheme 1 — Leave-One-Out (LOO-CV):**

Compounds are held out one at a time, with the model trained on the remaining
45. This tests whether the model can interpolate within the chemical space
of the dataset.

| Model | R² | MAE | Pearson r | p-value |
|---|---|---|---|---|
| Ridge | 0.708 | 0.136 | 0.843 | 1.92×10⁻¹³ |
| Lasso | 0.777 | 0.112 | 0.883 | 4.89×10⁻¹⁶ |
| Random Forest | 0.752 | 0.124 | 0.871 | 3.31×10⁻¹⁵ |
| Gradient Boosting | 0.763 | 0.121 | 0.874 | 2.10×10⁻¹⁵ |

**Cross-validation scheme 2 — Leave-One-Series-Out (LOSO-CV):**

Each chemical series is held out in turn while the model trains on the other
two. This tests whether the model generalises *across* chemical scaffolds —
i.e., whether it has learned genuine structure-activity relationships or
merely learned to identify which series a compound belongs to.

| Model | R² | MAE | Pearson r | p-value |
|---|---|---|---|---|
| Ridge | −0.182 | 0.319 | 0.074 | 0.627 |
| Lasso | −0.385 | 0.342 | 0.006 | 0.970 |
| Random Forest | −0.139 | 0.324 | −0.002 | 0.991 |
| Gradient Boosting | −0.377 | 0.342 | −0.027 | 0.857 |

LOSO R² is negative for all models, meaning they perform worse than
predicting the mean. The LOO-LOSO R² drop (~0.9 for linear models) confirms
that the strong LOO performance was entirely an artefact of inter-series
scaffold confounding: the model learned the identity of each chemical series
(and its characteristic bias value) rather than the within-series
structure-activity relationship.

**Cross-validation scheme 3 — Within-series LOO-CV:**

To test whether genuine SAR exists within each scaffold independently,
separate models were trained and evaluated within each series using LOO-CV.

| Series | n | Ridge R² | Ridge p-value | Signal? |
|---|---|---|---|---|
| Heterocyclic | 17 | 0.413 | 5.2×10⁻³ | Weak–moderate |
| Pottie 2023 | 17 | 0.452 | 9.4×10⁻⁴ | Weak–moderate |
| 25CN-NBOH | 12 | 0.328 | 2.7×10⁻² | See below |

**Sensitivity analysis — 25CN-NBOH series:** Serotonin is included in this
series as a reference compound from the source paper, but is structurally
anomalous: it is an indoleamine (not an N-benzyl phenethylamine) and has the
most negative log_bias in the series (−0.305). Removing it collapses R² from
0.328 to 0.024, revealing that the apparent signal was a leverage effect from
one structurally distinct outlier. The 25CN-NBOH series has no reliable
within-series SAR signal at n=11.

### Summary of findings

1. **LOO-CV is misleadingly optimistic** for datasets with few, chemically
   distinct series. The apparent R²=0.78 was inflated by scaffold confounding
   and should not be reported without the LOSO result alongside it.

2. **Weak but genuine within-series SAR signal exists** in the Heterocyclic
   and Pottie 2023 series (Ridge R²≈0.41–0.45, p<0.01), meaning Morgan
   fingerprints capture some structural features that correlate with bias
   within a scaffold. The signal is too weak for reliable predictions but
   is above chance.

3. **Cross-scaffold generalisation fails completely.** A model trained on two
   series cannot predict the third. This is not surprising: with only three
   chemical series and ~17 compounds each, there is insufficient data to
   learn generalizable structure-bias relationships from 2D fingerprints.

4. **The primary bottleneck is data, not methodology.** A meaningful
   cross-series QSAR model would require at minimum 5–10 distinct chemical
   scaffolds with paired Gq/β-arrestin measurements on the same assay
   platform, with at least 10–15 compounds per scaffold. ChEMBL currently
   has very limited β-arrestin assay coverage for 5-HT2A (only 7 assays
   total, from 3 papers).

### Reproducing the analysis

Dependencies:

```bash
conda activate molml  # or any environment with rdkit, scikit-learn, pandas
pip install requests pandas matplotlib scikit-learn
# rdkit should be installed via conda: conda install -c conda-forge rdkit
```

Run scripts in order:

```bash
python analysis/01_explore_chembl_functional.py   # fetch ChEMBL functional data
python analysis/02_identify_chembl_sources.py     # identify source papers
python analysis/03_compute_bias_factors.py        # compute log bias for ChEMBL compounds
python analysis/04_merge_and_pubchem.py           # merge with literature data
python analysis/05_qsar_model.py                  # LOO-CV and LOSO-CV
python analysis/06_within_series_qsar.py          # within-series models + sensitivity
```

Script 04 requires `analysis/literature_curated.csv`, which contains manually
transcribed data from Table 1 of Pottie & Poulie et al. (ACS Chem. Neurosci.
2023, DOI: 10.1021/acschemneuro.3c00267). The transcription is included in
this repository; the original paper is open access at
https://pmc.ncbi.nlm.nih.gov/articles/PMC10401645/.

Key output files:

| File | Description |
|---|---|
| `analysis/assay_summary.csv` | All 295 functional assays for 5-HT2A with pathway classification |
| `analysis/paired_compounds.csv` | 29 ChEMBL compounds with paired Gq/βarr Emax |
| `analysis/bias_factors.csv` | Log bias factors for ChEMBL compounds |
| `analysis/literature_curated.csv` | Manually curated data from Pottie & Poulie 2023 |
| `analysis/combined_dataset.csv` | Full 46-compound dataset with SMILES and log bias |
| `analysis/bias_combined.png` | Bias distribution across all sources |
| `analysis/qsar_loo_predictions.png` | LOO-CV predicted vs actual |
| `analysis/qsar_loso_predictions.png` | LOSO-CV predicted vs actual |
| `analysis/within_series_predictions.png` | Within-series LOO-CV plots |
| `analysis/within_series_sensitivity.png` | Sensitivity analysis (serotonin removed) |
| `analysis/05_qsar_report.txt` | Full LOO and LOSO results |
| `analysis/06_within_series_report.txt` | Within-series results |
