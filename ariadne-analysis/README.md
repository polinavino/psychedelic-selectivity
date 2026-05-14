# Ariadne Pharmacological Analysis

Quantitative re-analysis of data from Cunningham et al. ACS Chem Neurosci
2023 (doi:10.1021/acschemneuro.2c00597), applying bias factor and selectivity
quantification to the published numerical values.

## Background

Ariadne (4C-D, BL-3912A) is a phenylalkylamine synthesised by Alexander
Shulgin, structurally identical to the hallucinogen DOM except for a single
additional methylene group at the alpha-position (alpha-ethyl vs alpha-methyl).
Despite being a potent 5-HT2A agonist, Ariadne produces no hallucinogenic
effects in humans at doses up to 300 mg, while showing remarkable therapeutic
effects in Bristol-Myers clinical trials in the 1970s: rapid remission of
psychotic symptoms, near-complete remission of Parkinson's symptoms, and
pro-cognitive effects in geriatric subjects.

Cunningham et al. 2023 characterised its pharmacology and proposed the
**signalling efficacy hypothesis**: Ariadne is non-hallucinogenic not because
of biased agonism but because the alpha-ethyl group causes steric clash with
residues W336/F340 in the orthosteric binding pocket of 5-HT2A, uniformly
reducing potency and efficacy across all signalling pathways (Gq, G11,
Ca²⁺ flux, beta-arrestin2).

## What this analysis adds

The source paper is experimental; this analysis applies formal quantitative
frameworks to the published numbers. Three analyses are included.

### Analysis 1 — Assay-dependent efficacy (ariadne_gq_vs_caflux.png)

Gq dissociation BRET and Ca²⁺ flux results are presented side-by-side in a
2×3 grid (rows = assay, columns = 5-HT2A/2B/2C), with pEC50 bars and Emax
annotated on a twin axis.

**Key finding:** At 5-HT2C, Ariadne's Emax drops from ~82% (Gq BRET) to
~34% (Ca²⁺ flux) for the racemate, and from ~79% to ~22% for (S)-Ariadne —
a 2-4 fold difference in apparent efficacy from assay choice alone. Potency
(EC50) is similar across assays; only efficacy diverges, and only at 5-HT2C.
This visualisation is not present in the original paper but is directly
computable from its Figure 2d data.

The practical implication: the apparent selectivity profile of Ariadne
across 5-HT2 subtypes depends heavily on which assay is used. A researcher
relying on Gq BRET would conclude Ariadne is a moderate 5-HT2C agonist
(Emax ~82%); Ca²⁺ flux gives a very different picture (Emax ~34%).

### Analysis 2 — Bias factor across alpha-alkyl pairs (ariadne_alpha_series.png)

The log bias factor is computed for three matched hallucinogenic/non-hallucinogenic
pairs at 5-HT2A:

log_bias = log10( (Emax_barr / EC50_barr) / (Emax_Gq / EC50_Gq) )

| Compound | Alpha | log_bias | Hallucinogenic? |
|---|---|---|---|
| DOM | methyl | −0.642 | YES |
| (rac)-Ariadne | ethyl | −0.542 | no |
| DOPr | methyl | −0.721 | YES |
| 4C-Pr | ethyl | −0.368 | no |
| DOI | methyl | −0.735 | YES |
| 4C-I | ethyl | −0.449 | no |

**Key finding:** All six compounds are Gq-biased (negative log_bias). The
alpha-ethyl (non-hallucinogenic) compound in each pair is consistently
slightly less Gq-biased than its alpha-methyl (hallucinogenic) counterpart,
by 0.10-0.35 log units. The direction is consistent across all three
4-position substitutions but the magnitude is small.

This formally confirms the paper's signalling efficacy hypothesis: Ariadne
is not a biased agonist. Both pathways drop together. The mechanism of
non-hallucinogenicity is uniform partial agonism, not pathway selectivity.

The paper discusses each pair in a separate figure panel; this is the first
unified quantitative comparison across all three pairs.

### Analysis 3 — 4-position analog 2A/2B selectivity (ariadne_2a2b_selectivity.png)

pEC50 differences (pEC50_2A - pEC50_2B) computed for all five 4-position
analogs reported in Figure 5 of the paper:

| Compound | pEC50 2A-2B | Fold selectivity |
|---|---|---|
| 4C-D (Ariadne) | +0.67 | 4.7x |
| 4C-TFM | +0.91 | 8.2x |
| 4C-Pr | +0.91 | 8.2x |
| 4C-cycPr | +0.69 | 4.9x |
| 4C-MOM | −0.02 | ~1x (balanced) |

**Key finding:** 4C-TFM and 4C-Pr show the highest 5-HT2A/2B selectivity
(+0.91, 8-fold each). 4C-MOM is essentially balanced and has substantially
reduced Emax at 5-HT2B (43%), likely due to increased polarity. The paper
identifies 4C-TFM as most selective but does not express selectivity as
pEC50 differences or compare directly across all five analogs in one table.

2A/2B selectivity is clinically relevant: sustained 5-HT2B agonism causes
cardiac valvulopathy (the fenfluramine mechanism). All Ariadne analogs
except 4C-MOM are 2A-preferring, which is favourable for therapeutic use.

## What is not included and why

**Four-metric selectivity comparison (S-score, Gini, entropy, ratio)** is
not computed. With only 3 receptors (5-HT2A/2B/2C) having tabular EC50/Emax
data in the paper, all selectivity definitions are mathematically forced to
agree — there is no statistical room for divergence. The entropy vs Gini
finding that replicates across kinase and PDSP datasets requires broader
receptor coverage to emerge. This analysis would require Table S3 from the
paper (full 12-receptor 5-HTome numerical data), which is referenced in the
paper but not in the publicly available PMC version.

## Data

All values are transcribed from Cunningham et al. 2023 figures and tables:

- Figure 2d: 5-HT2A/2B/2C Gq BRET and Ca²⁺ flux EC50/Emax for 5-HT and
  Ariadne enantiomers
- Figure 2a: Ki values (radioligand displacement) at 5-HT2A
- Figures 3a/b/c: DOM vs Ariadne, DOPr vs 4C-Pr, DOI vs 4C-I at 5-HT2A
  (Gq and beta-arrestin2)
- Figure 3d: DOM vs Ariadne vs 5C-D at 5-HT2A (Gq only; no beta-arrestin2
  data for 5C-D in the paper)
- Figure 5c: 4-position analogs at 5-HT2A/2B/2C (Gq)

See `data/cunningham2023_data_notes.txt` for column definitions, assay
details, and compound nomenclature.

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
| `PRIMARY_ASSAY` | gq_bret | Assay for potency profiles |

### Outputs

| File | Description |
|---|---|
| `analysis/ariadne_report.txt` | Full numerical results |
| `analysis/ariadne_gq_vs_caflux.png` | Assay-dependent efficacy comparison |
| `analysis/ariadne_alpha_series.png` | Bias factor across alpha-alkyl pairs |
| `analysis/ariadne_2a2b_selectivity.png` | 4-position analog 2A/2B selectivity |

## Connection to companion work

- `../analysis/` — biased agonism QSAR at 5-HT2A with LOSO evaluation
- `../shulgin-selectivity/` — PDSP binding profiles for 36 Shulgin compounds
  including DOM and DOI (the hallucinogenic counterparts of Ariadne and 4C-I)
- `../../kinase-selectivity-definitions/` — the selectivity definition
  framework applied here
