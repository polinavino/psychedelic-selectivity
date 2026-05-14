# Ariadne Selectivity Analysis

Quantitative pharmacological analysis of Ariadne (4C-D, BL-3912A) and its
structural analogs, applying the serotonin receptor selectivity framework
from the companion analysis to data from Cunningham et al. 2023.

## Background

Ariadne is a phenylalkylamine synthesised by Alexander Shulgin, structurally
identical to the hallucinogen DOM except for a single additional methylene
group at the alpha-position (alpha-ethyl vs alpha-methyl). Despite being a
potent 5-HT2A agonist, Ariadne produces no hallucinogenic effects in humans
at doses up to 300 mg, while showing remarkable therapeutic effects including
rapid remission of psychotic symptoms, anti-Parkinsonian effects, and
pro-cognitive effects in clinical trials conducted by Bristol-Myers in the
1970s. The compound was abandoned and its molecular pharmacology remained
unknown for decades.

Cunningham et al. (ACS Chem Neurosci 2023, doi:10.1021/acschemneuro.2c00597)
characterised Ariadne's pharmacology in detail and proposed the **signalling
efficacy hypothesis**: Ariadne is non-hallucinogenic not because of biased
agonism (it is not), but because it is a weaker partial agonist at 5-HT2A
across all signalling pathways simultaneously (Gq, G11, Ca²⁺, beta-arrestin2).
The alpha-ethyl group causes steric clash with W336/F340 in the orthosteric
binding sub-pocket, reducing potency and efficacy for everything.

## What this analysis adds

The Cunningham 2023 paper reports the experimental pharmacology and a
qualitative docking rationale. This analysis applies formal selectivity
and bias quantification frameworks — developed in the companion kinase
selectivity and serotonin receptor analyses — to the published numerical data.

Specifically:

**1. Selectivity metrics across 5-HT2A/2B/2C**
S-score, Gini coefficient, Shannon entropy, and selectivity ratio are computed
for Ariadne and its enantiomers and compared against 5-HT as reference. This
is the first application of this four-metric framework to the Ariadne series.

**2. Assay-dependent efficacy profile (Gq BRET vs Ca²⁺ flux)**
The two assays yield dramatically different Emax values at 5-HT2C: Ariadne
shows ~82% Emax in Gq BRET but only ~34% in Ca²⁺ flux. This assay dependence
is visualised side-by-side in a way the original paper does not present,
highlighting how assay choice affects apparent receptor selectivity.

**3. Bias factor quantification across the alpha-alkyl series**
The log_bias = log10((Emax_barr/EC50_barr) / (Emax_Gq/EC50_Gq)) is computed
for each matched hallucinogenic/non-hallucinogenic pair (DOM/Ariadne,
DOPr/4C-Pr, DOI/4C-I). The original paper describes the shift qualitatively;
this analysis quantifies it across all three pairs simultaneously.

**4. 4-position analog 2A/2B selectivity**
pEC50 differences (pEC50_2A - pEC50_2B) are computed for all 4-position
analogs of Ariadne reported in Figure 5 of the paper. 4C-TFM and 4C-Pr
show the highest 2A/2B selectivity (+0.91, corresponding to 8-fold).

## Key findings

### Selectivity metrics

With only 3 receptors measured, all four selectivity definitions agree
(r > 0.96 in all pairwise correlations). This contrasts with the kinase
and PDSP serotonin receptor datasets where entropy diverges from Gini —
consistent with the theoretical prediction that the divergence requires
sufficient receptor coverage to reveal distributional differences. With
n=3 receptors, the distributions are too constrained.

5-HT has higher Gini than Ariadne (0.040 vs 0.020-0.024), reflecting its
dramatically higher 5-HT2C potency relative to 2A/2B. Ariadne's Gini values
are lower because its potency is more uniformly distributed across the three
5-HT2 subtypes — an observation absent from the original paper.

### Bias factors

| Compound | Alpha | log_bias | Hallucinogenic? |
|---|---|---|---|
| DOM | methyl | −0.642 | YES |
| (rac)-Ariadne | ethyl | −0.542 | no |
| DOPr | methyl | −0.721 | YES |
| 4C-Pr | ethyl | −0.368 | no |
| DOI | methyl | −0.735 | YES |
| 4C-I | ethyl | −0.449 | no |

All compounds are Gq-biased (negative log_bias). The alpha-ethyl (Ariadne)
series consistently shows slightly less Gq bias than the alpha-methyl (DOx)
hallucinogens, by approximately 0.1-0.35 log units. The direction is
consistent across all three 4-position substitution pairs but the magnitude
is small — confirming the paper's conclusion that the mechanism is not
biased agonism but rather a uniform reduction in efficacy across all pathways.

### 4-position analog selectivity

| Compound | pEC50 2A-2B | Notes |
|---|---|---|
| 4C-D (Ariadne) | +0.67 | baseline |
| 4C-TFM | +0.91 | 8-fold selectivity — most selective |
| 4C-Pr | +0.91 | equally selective, lower Emax at 2B |
| 4C-cycPr | +0.69 | similar to baseline |
| 4C-MOM | −0.02 | essentially balanced; highest 2B Emax loss |

4C-TFM is the most therapeutically promising analog: high 2A/2B selectivity
(important for cardiac safety), non-hallucinogenic (alpha-ethyl), and potent
at 5-HT2A (EC50 = 28.6 nM). The paper identifies it but does not express its
selectivity as a pEC50 difference or compare it to the other analogs in this
format.

## Data

All data are transcribed from figures and tables in Cunningham et al. 2023.
The data file includes:

- Figure 2d: 5-HT2A/2B/2C Gq dissociation BRET and Ca²⁺ flux EC50/Emax
  for 5-HT, (rac)-, (R)-, and (S)-Ariadne
- Figure 3a: DOM vs Ariadne at 5-HT2A, Gq and beta-arrestin2
- Figure 3b: DOPr vs 4C-Pr at 5-HT2A, Gq and beta-arrestin2
- Figure 3c: DOI vs 4C-I at 5-HT2A, Gq and beta-arrestin2
- Figure 3d: DOM vs Ariadne vs 5C-D at 5-HT2A, Gq
- Figure 5c: 4-position analogs at 5-HT2A/2B/2C, Gq
- Figure 2a: Ki values (radioligand displacement) for Ariadne enantiomers

See `data/cunningham2023_data_notes.txt` for full column definitions, assay
details, and nomenclature.

## Reproducing the analysis

```bash
cd psychedelic-selectivity/ariadne-analysis/
python analysis/ariadne_selectivity.py
```

### Dependencies

```bash
conda activate molml
pip install pandas matplotlib numpy
```

### Parameters

| Parameter | Default | Description |
|---|---|---|
| `SSCORE_THRESHOLD` | 6.5 | pEC50 threshold for S-score (~316 nM) |
| `PRIMARY_ASSAY` | gq_bret | Assay used for selectivity profiles |

### Outputs

| File | Description |
|---|---|
| `analysis/ariadne_selectivity_report.txt` | Full numerical results |
| `analysis/ariadne_gq_vs_caflux.png` | 2x3 grid: Gq BRET vs Ca²⁺ flux potency and efficacy |
| `analysis/ariadne_alpha_series.png` | Alpha-alkyl series bias across three 4-position pairs |
| `analysis/ariadne_profiles.png` | Overview: pEC50 profile and alpha-chain trend |
| `analysis/ariadne_selectivity_scores.png` | Four selectivity metrics comparison |
| `data/ariadne_series.csv` | Computed pEC50 and log_transduction values |

## Limitations

- Only 3 receptors (5-HT2A, 5-HT2B, 5-HT2C) have tabular EC50/Emax data
  in the paper. The 5-HTome screen (Figure 2b) covers all 12 serotonin
  receptors but reports only log(Emax/EC50) as a colour scale, not numerical
  values. A full 12-receptor selectivity profile would require contacting the
  authors for supplementary data (Table S3 mentioned in the paper but not
  publicly available as of this writing).
- With only 3 receptors, the selectivity metric correlations are not
  meaningful — all definitions agree trivially. The entropy divergence finding
  from the kinase and PDSP datasets requires more receptors to emerge.
- Beta-arrestin2 data is available for DOM and Ariadne (Fig 3a) and two
  matched pairs (Figs 3b, 3c) but not for 5C-D — so log_bias cannot be
  computed for the alpha-propyl compound.
- All data are from HEK293T cells transiently transfected with human
  receptors. Native tissue pharmacology may differ.
- Enantiomers: (R)-Ariadne is the therapeutically active enantiomer
  (used in Bristol-Myers trials). The racemate is approximately the average
  of R and S in most assays.

## Connection to companion work

This analysis is part of a larger project applying formal selectivity
definitions to serotonin receptor pharmacology:

- `../analysis/` — biased agonism QSAR at 5-HT2A (LOO vs LOSO evaluation)
- `../shulgin-selectivity/` — PDSP Ki database selectivity profiles for 36
  Shulgin compounds including DOM, DOI, and related compounds
- `../../kinase-selectivity-definitions/` — the original four-metric
  selectivity framework (pre-print: https://chemrxiv.org/doi/full/10.26434/chemrxiv.15001618/v1)

The finding that Ariadne's Gini selectivity (0.020-0.024) is lower than
5-HT's (0.040) across 5-HT2A/2B/2C is consistent with the PDSP analysis
showing that most Shulgin phenethylamines bind the three 5-HT2 subtypes
with similar affinity — the alpha-chain extension in Ariadne uniformly
reduces this already-broad binding profile rather than sharpening it.
