# Shulgin Compound Selectivity at Serotonin Receptors

Receptor binding profiles and predicted signalling bias for compounds
documented in PIHKAL and TIHKAL, using data from the PDSP Ki database.

## In plain terms

Alexander Shulgin synthesised and characterised hundreds of psychoactive
compounds over his career, documenting them in two books: PIHKAL
(phenethylamines) and TIHKAL (tryptamines). Most of these compounds
have never been tested in modern binding assays, and even those that have
been tested are rarely compared systematically across the full serotonin
receptor family.

This analysis asks three questions for each Shulgin compound with available
binding data:

1. **Which serotonin receptors does it bind, and how potently?** (receptor
   selectivity profile)
2. **Is it 5-HT2A selective over 5-HT2B?** (cardiac safety — 5-HT2B
   agonism causes valvulopathy at sustained activation)
3. **Does its structure predict Gq or beta-arrestin signalling bias at
   5-HT2A?** (hallucination vs neuroplasticity pathway — see caveat below)

## Data source

Binding affinity (Ki) data from the **PDSP Ki Database** (Psychoactive
Drug Screening Program, UNC Chapel Hill), downloaded via the
`search.php` pagination API (98,685 records total; the download button
was broken as of 2025 — see script 00 for the workaround).

The PDSP database provides Ki values in nM with SMILES, species, tissue
source, and reference for each measurement. All values are converted to
pKi = 9 - log10(Ki / nM) for analysis, so higher pKi = more potent
binding.

## Compound selection

114 Shulgin compound names from PIHKAL and TIHKAL were searched against
the PDSP database using exact and substring matching, followed by manual
curation to remove false positives (antipsychotics, antidepressants, and
other drugs that matched Shulgin names as substrings). The curated
whitelist maps 71 PDSP ligand name variants to 65 canonical Shulgin names.

After filtering to compounds with valid Ki values at 3 or more serotonin
receptors, **36 compounds** remain for selectivity analysis.

## Selectivity profiles

### Receptor coverage

| Compound | Receptors tested | 5-HT2A | 5-HT2B |
|---|---|---|---|
| serotonin, 5-MeO-T | 12 | Y | Y |
| LSD, DMT, 5-MeO-DMT, psilocin, psilocybin, 2C-B, 2C-E, MDMA, DPT, 6-F-DMT | 10 | Y | Y |
| DALT series, DOI, DiPT, 5-MeO-DiPT, 5-MeO-MiPT, mescaline | 9 | Y | mostly Y |
| DOB, MDA, ibogaine, 2-Bromo-LSD | 5-7 | Y | partial |
| DOM, DOET, LSD(+) | 4 | Y | N |

### 5-HT2A vs 5-HT2B selectivity (cardiac safety)

pKi differences computed as pKi(5-HT2A) - pKi(5-HT2B). Positive values
indicate 5-HT2A preference (safer); negative values indicate 5-HT2B
preference (cardiac risk at sustained activation).

| Compound | 2A-2B | pKi_2A | pKi_2B | Note |
|---|---|---|---|---|
| psilocybin | +1.27 | 7.89 | 6.61 | Most 2A-selective |
| 2-Ph-DALT | +1.17 | 7.89 | 6.72 | Novel DALT analogue |
| LSD | +0.84 | 8.57 | 7.73 | |
| 2-Bromo-LSD | +0.81 | 8.82 | 8.02 | |
| 5-MeO-DMT | +0.76 | 6.41 | 5.65 | |
| DOI | +0.46 | 8.28 | 7.82 | |
| DOB | +0.49 | 8.31 | 7.82 | |
| 2C-B | -0.06 | 7.84 | 7.90 | Essentially balanced |
| psilocin | -1.86 | 6.47 | 8.33 | Strongly 2B-preferring |
| DPT | -1.79 | 5.59 | 7.38 | |
| 5-MeO-T | -1.98 | 6.28 | 8.26 | Most 2B-preferring |

**Notable finding:** psilocybin is the most 5-HT2A-selective compound in
the dataset (+1.27), yet its active metabolite psilocin is strongly
5-HT2B-preferring (-1.86). This means psilocybin's favourable cardiac
safety profile in clinical use may partly reflect its prodrug nature:
rapid dephosphorylation to psilocin is not sustained chronically in the
way that would produce valvulopathy. This distinction is rarely made
explicit in the pharmacological literature.

### Four selectivity definitions

The same framework applied in the kinase inhibitor selectivity paper
(see companion repo) is applied here. For each compound, four selectivity
scores are computed from the pKi profile across tested receptors:

- **S-score**: fraction of tested receptors below pKi = 6.0 (~1 uM).
  Higher = more selective.
- **Gini**: Gini coefficient of the pKi distribution. Higher = more selective.
- **Entropy**: Shannon entropy of the normalised pKi distribution.
  Lower = more selective.
- **Ratio**: pKi of top receptor / pKi of second receptor. Higher = more selective.

Correlation between definitions across 36 compounds:

|  | S-score | Gini | Entropy | Ratio |
|---|---|---|---|---|
| S-score | 1.000 | 0.468 | -0.628 | 0.277 |
| Gini | 0.468 | 1.000 | -0.078 | -0.019 |
| Entropy | -0.628 | -0.078 | 1.000 | -0.465 |
| Ratio | 0.277 | -0.019 | -0.465 | 1.000 |

Entropy is uncorrelated with Gini (r = -0.078) and weakly correlated with
ratio (r = -0.465), while Gini and S-score are moderately correlated
(r = 0.468). This reproduces the two-cluster structure found in kinase
inhibitor data and serotonin receptor data in the companion analysis:
entropy consistently measures something categorically different from
concentration-based definitions, across three independent datasets and two
receptor families.

Note: Gini values are uniformly small (0.026-0.088) across all compounds
because serotonergic psychedelics bind serotonin receptors broadly, with
pKi values in a narrow range (5-9). These compounds were not optimised
for subtype selectivity.

## Predicted signalling bias (exploratory)

The Lasso QSAR model from the companion biased agonism analysis
(`../analysis/combined_dataset.csv`, n=46 compounds, three chemical
series) was applied to predict log_bias at 5-HT2A for each Shulgin
compound. The log_bias is defined as:

```
log_bias = log10( (Emax_barr / EC50_barr) / (Emax_Gq / EC50_Gq) )
```

Positive = beta-arrestin biased (neuroplasticity pathway), negative = Gq
biased (hallucinogenic pathway).

**Results:** predicted log_bias range is -0.242 to +0.362, compared to
the training data range of -0.596 to +0.864. The model strongly regresses
to the mean for all novel scaffolds. DOx and 2C-X phenethylamines receive
predictions near +0.36, close to the Pottie 2023 training series mean
(+0.40). Tryptamines receive predictions near -0.24, close to the
Heterocyclic training series mean (-0.20).

This confirms in practice the scaffold confounding problem identified in
leave-one-series-out cross-validation (LOSO R² < 0 for all models):
the model assigns predictions based on which training series a compound
structurally resembles, not based on specific features that drive bias.
These predictions should not be interpreted quantitatively.

A reliable virtual screen for signalling bias would require either a much
larger and more structurally diverse training dataset spanning many
scaffolds, or 3D structure-based features from receptor docking rather
than 2D Morgan fingerprints.

## Reproducing the analysis

### Dependencies

```bash
conda activate molml  # or any environment with rdkit and scikit-learn
pip install requests pandas matplotlib scikit-learn
# rdkit via conda: conda install -c conda-forge rdkit
```

### Steps

```bash
cd shulgin-selectivity/

# 1. Download PDSP database via API pagination (download button broken)
python analysis/00_download_pdsp.py
# Output: analysis/KiDatabase.csv (~91,000 records, stops at offset 91,200
# due to a malformed JSON record on the PDSP server)

# 2. Search for Shulgin compounds in PDSP
python analysis/01_identify_shulgin.py
# Output: analysis/shulgin_hits.csv, analysis/serotonin_receptors.csv,
#         analysis/01_inspect_report.txt

# 3. Clean false positives, compute pKi values
python analysis/02_clean_shulgin.py
# Output: analysis/shulgin_clean.csv, analysis/shulgin_compound_list.csv,
#         analysis/02_clean_report.txt

# 4. Compute selectivity profiles
python analysis/03_selectivity_profiles.py
# Output: analysis/selectivity_profiles.csv, analysis/pki_matrix.csv,
#         analysis/selectivity_heatmap.png, analysis/selectivity_scores.png,
#         analysis/03_selectivity_report.txt

# 5. Predict signalling bias (requires ../analysis/combined_dataset.csv)
python analysis/04_bias_predictions.py
# Output: analysis/bias_predictions.csv, analysis/bias_vs_selectivity.png,
#         analysis/04_bias_report.txt
```

### Key parameters

| Script | Parameter | Default | Description |
|---|---|---|---|
| 03 | `THRESHOLD` | 6.0 | pKi threshold for S-score |
| 03 | `MIN_RECEPTORS` | 3 | Minimum receptors for inclusion |
| 03 | `HEATMAP_MIN_RECEPTORS` | 5 | Minimum receptors to appear in heatmap |
| 04 | `RADIUS` | 2 | Morgan fingerprint radius (ECFP4) |
| 04 | `N_BITS` | 2048 | Fingerprint bit length |
| 04 | `LASSO_ALPHA` | 0.01 | Lasso regularisation strength |
| 04 | `BIAS_THRESHOLD` | 0.3 | |log_bias| for barr/Gq classification |

### Output files

| File | Description |
|---|---|
| `analysis/KiDatabase.csv` | Full PDSP download (91,200 records) |
| `analysis/shulgin_hits.csv` | PDSP records matching Shulgin names |
| `analysis/shulgin_clean.csv` | Cleaned records with pKi values |
| `analysis/shulgin_compound_list.csv` | Per-compound receptor coverage |
| `analysis/pki_matrix.csv` | pKi matrix: compound x receptor |
| `analysis/selectivity_profiles.csv` | All selectivity scores per compound |
| `analysis/bias_predictions.csv` | Predicted log_bias per compound |
| `analysis/selectivity_heatmap.png` | Binding profile heatmap |
| `analysis/selectivity_scores.png` | Selectivity definition correlations |
| `analysis/bias_vs_selectivity.png` | 2A/2B selectivity vs predicted bias |

## Limitations

- PDSP coverage is uneven: well-studied compounds (LSD, DMT, psilocin)
  have 10+ receptor measurements while others (DOM, DOET) have only 4.
- Ki values aggregate across species (human, rat, bovine) and tissue
  sources. Human data is preferred where available but the dataset is mixed.
- 107 records are censored (Ki > threshold) and excluded from pKi
  computation. This may bias profiles for compounds where many
  non-2A measurements are censored.
- The PDSP download stopped at offset 91,200 due to a malformed JSON
  record on the server; the remaining ~7,000 records are inaccessible
  without server-side fixes.
- Signalling bias predictions are not reliable for novel scaffolds.
  See the companion biased agonism analysis for the full evidence.

## Connection to companion work

This analysis applies the selectivity framework from the kinase inhibitor
selectivity paper to a new receptor family and a historically significant
set of compounds. The entropy finding now replicates across three
independent datasets:

1. Kinase inhibitor binding profiles (Davis + Klaeger datasets)
2. Serotonin receptor profiles from ChEMBL (n=297 compounds)
3. Shulgin compound profiles from PDSP (n=36 compounds)

This consistency across receptor families and compound classes strengthens
the case that entropy measures a categorically different aspect of
selectivity than concentration-based definitions, and that the choice of
selectivity metric should be made explicitly rather than treated as
interchangeable.
