# 5-HT2A Molecular Docking — Shulgin Compounds

AutoDock Vina docking of 41 Shulgin compounds into the active-state
5-HT2A crystal structure (PDB: 6WHA), with correlation analysis against
measured binding affinities from the PDSP Ki database.

## In plain terms

Molecular docking predicts how a small molecule binds to a protein by
computationally searching for the lowest-energy binding pose. A good
docking model should give more negative scores (better predicted binding)
for compounds that are known to bind tightly, and less negative scores for
weak binders.

This analysis tests whether AutoDock Vina docking scores at the 6WHA
structure correlate with measured 5-HT2A pKi values from the PDSP database.
The answer is no — and the reason is chemically interpretable.

## Structure

**PDB: 6WHA** — cryo-EM structure of 5-HT2A bound to the agonist 25CN-NBOH
in complex with mini-Gαq, Gβ1, Gγ2, and scFv16 (active state, 2020).

25CN-NBOH is a tryptamine-derived compound. Its binding induces a specific
conformation of the orthosteric pocket, particularly displacing W336 and
F340 to accommodate its phenoxy group. This open conformation of the
sub-pocket is not seen in 5-HT2B or 5-HT2C structures with non-tryptamine
ligands, where the pocket is partially collapsed.

**Chain extraction:** chain A residues 66-404 only (receptor). Chains B-E
(Gαq, Gβ1, Gγ2, scFv16) and the co-crystallised ligand U0G (25CN-NBOH)
are removed. Binding box centred on the U0G centroid (98.13, 106.21, 94.14),
25 Å cube.

## Results

### Docking scores

41/43 compounds docked successfully. LSD and 2-Bromo-LSD failed with a
Vina internal error — likely due to the ergoline scaffold's PDBQT rotatable
bond definition. All other compounds completed in 2-6 seconds at
exhaustiveness=8.

Score range: −8.66 (ibogaine) to −5.18 (2C-B) kcal/mol.

### Correlation with measured pKi

| Compound set | n | Pearson r | p-value | Direction |
|---|---|---|---|---|
| All | 31 | +0.428 | 0.016 | WRONG |
| Tryptamines | 22 | −0.243 | 0.277 | correct but n.s. |
| Phenethylamines | 9 | +0.800 | 0.010 | WRONG |

**Overall r = +0.428** is positive — the wrong direction. A predictive
docking model should show negative correlation (better docking score →
higher pKi). Instead, compounds that dock better (more negative score)
tend to have *lower* measured affinity.

### Scaffold-dependent docking performance (3D scaffold confounding)

The apparent positive correlation is driven entirely by scaffold separation:

| Scaffold class | Mean Vina score | Mean pKi at 5-HT2A |
|---|---|---|
| Tryptamines (n=22) | −6.67 kcal/mol | 6.34 |
| Phenethylamines (n=9) | −5.69 kcal/mol | 7.49 |

Tryptamines dock **0.97 kcal/mol better** than phenethylamines on average
(t=−6.07, p=1.3×10⁻⁶), yet phenethylamines have substantially higher measured
binding affinity. DOB and DOI have pKi = 8.27 and 8.28 — the strongest binders
in the dataset — but dock at only −5.57 and −5.66 kcal/mol. 2C-B (pKi=7.81)
is the weakest docking compound at −5.18 kcal/mol.

The within-phenethylamine correlation of r=+0.80 in the wrong direction
further illustrates this: even within the phenethylamine series, the model
ranks compounds in reverse order of true affinity.

**Cause:** 6WHA was co-crystallised with 25CN-NBOH, a tryptamine-like
compound. The orthosteric pocket geometry in this structure is optimised for
the tryptamine binding mode. Phenethylamines (2C-X, DOx series) bind to
5-HT2A through a different binding mode — they sit higher in the pocket and
make different contacts — and this mode is not captured by the 6WHA
conformation.

**Practical implication:** Virtual screening at 6WHA will systematically
rank tryptamines above phenethylamines regardless of true affinity. A
phenethylamine hit from a 6WHA-based screen would need to score below −6.5
kcal/mol to be competitive with tryptamines, when phenethylamines in this
dataset average only −5.7.

### Connection to the QSAR scaffold confounding paper

This is the 3D analogue of the scaffold confounding finding from the
companion ADMET QSAR analysis:

- 2D QSAR: Morgan fingerprint models trained on one scaffold (e.g.
  tryptamines) cannot predict activity for another scaffold
  (phenethylamines) — LOSO R²<0.
- 3D docking: a structure co-crystallised with a tryptamine-like ligand
  cannot rank phenethylamines correctly — overall r in the wrong direction.

In both cases, the method has learned the training scaffold rather than
generalizable structure-property or structure-activity relationships.

### LSD docking failure

LSD and 2-Bromo-LSD failed with `internal_error` in Vina. The ergoline
scaffold is a fused tetracyclic system — larger and more rigid than any
other compound in the set. The PDBQT conversion by OpenBabel may define
rotatable bonds in a way that causes Vina's internal search to fail.
This could be fixed by manually editing the PDBQT or using a different
ligand preparation tool (e.g. Meeko). Not addressed here.

## Reproducing the analysis

### Dependencies

```bash
conda activate molml
pip install vina
brew install open-babel   # or: conda install -c conda-forge openbabel
pip install pandas matplotlib scipy numpy
```

### Steps

```bash
cd psychedelic-selectivity/docking/

# 1. Prepare receptor (extract chain A, convert to PDBQT, compute box)
python analysis/01_prepare_receptor.py

# 2. Prepare ligands (SMILES -> 3D PDBQT via OpenBabel)
python analysis/02_prepare_ligands.py

# 3. Run docking (each ligand in isolated subprocess)
python analysis/03_run_docking.py

# 4. Correlate with PDSP pKi values
python analysis/04_analyze_results.py

# 5. Summary figure
python analysis/05_summary_figure.py
```

### Key parameters

| Script | Parameter | Default | Description |
|---|---|---|---|
| 01 | `BOX_SIZE` | 25.0 | Docking box size in Angstrom |
| 03 | `EXHAUSTIVENESS` | 8 | Vina search thoroughness |
| 03 | `N_POSES` | 9 | Poses per ligand |
| 03 | `ENERGY_RANGE` | 3.0 | Max energy range from best pose |

### Output files

| File | Description |
|---|---|
| `data/6WHA.pdb` | Raw PDB download |
| `data/6WHA_receptor.pdb` | Chain A residues 66-404, no ligand |
| `data/6WHA_receptor.pdbqt` | Vina-ready receptor |
| `data/binding_site.txt` | Docking box parameters |
| `data/ligands/*.pdbqt` | 3D ligand files (45 compounds) |
| `data/ligands/ligand_manifest.csv` | Ligand preparation status |
| `results/docking_scores.csv` | Best Vina score per compound |
| `results/poses/*.pdbqt` | All docking poses per compound |
| `results/docking_pki_merged.csv` | Merged scores + pKi values |
| `analysis/docking_summary.png` | Three-panel summary figure |
| `analysis/docking_vs_pki.png` | Detailed scatter + box plots |

## Limitations

- Only one crystal structure used (6WHA). Different structures may give
  different results — in particular, a phenethylamine-bound structure would
  likely reverse the scaffold bias.
- Exhaustiveness=8 is Vina's default. Higher values (16, 32) give more
  reliable poses but were not needed for this correlation analysis.
- Protonation states assigned at pH 7.4 by OpenBabel. Actual binding-site
  protonation may differ.
- Docking does not account for receptor flexibility. The binding pocket
  conformation is fixed at the 6WHA state (25CN-NBOH-bound active state).
- LSD and 2-Bromo-LSD could not be docked due to a Vina internal error
  with the ergoline scaffold.

## Connection to companion work

This analysis is part of a series:

- `../analysis/` — biased agonism QSAR at 5-HT2A (LOO vs LOSO scaffold
  confounding in 2D)
- `../shulgin-selectivity/` — PDSP pKi values used here for correlation
- `../../admet-qsar-evaluation/` — LOSO scaffold confounding in ADMET QSAR
  across ESOL, Lipophilicity, hERG (2D)
- `../../kinase-selectivity-definitions/` — selectivity definition framework
