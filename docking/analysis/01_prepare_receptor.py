"""
Script 01: Prepare 5-HT2A receptor from PDB 6WHA for AutoDock Vina.

6WHA is a cryo-EM structure of 5-HT2A (HTR2A) bound to the agonist
25CN-NBOH (U0G) in complex with mini-Gaq, Gb1, Gg2, and scFv16.

Chain map:
  Chain A residues -42 to 62 : E.coli C562 fusion (remove)
  Chain A residues 66 to 404 : 5-HT2A receptor (KEEP)
  Chain B                    : mini-Galpha-q (remove)
  Chain C                    : Gbeta1 (remove)
  Chain D                    : Ggamma2 (remove)
  Chain E                    : scFv16 antibody fragment (remove)
  HETATM U0G                 : 25CN-NBOH co-crystallised agonist (remove,
                               but use its centroid for docking box)

Steps:
  1. Extract chain A residues 66-404
  2. Remove HETATM records
  3. Compute binding site centroid from U0G coordinates
  4. Convert to PDBQT with OpenBabel (pH 7.4 hydrogens, rigid receptor)

Outputs:
  data/6WHA_receptor.pdb      - cleaned receptor PDB
  data/6WHA_receptor.pdbqt    - Vina-ready receptor
  data/binding_site.txt       - docking box centre and size

Usage:
  cd psychedelic-selectivity/docking/
  python analysis/01_prepare_receptor.py
"""

import subprocess
import os

RAW_PDB        = "data/6WHA.pdb"
RECEPTOR_PDB   = "data/6WHA_receptor.pdb"
RECEPTOR_PDBQT = "data/6WHA_receptor.pdbqt"
LIGAND_RESNAME = "U0G"
BOX_SIZE       = 25.0   # Angstrom — large enough for full orthosteric pocket

# ── 1. Extract receptor chain A residues 66-404 ───────────────────────────────

print("Extracting chain A residues 66-404...")

receptor_lines = []
ligand_lines   = []

with open(RAW_PDB) as f:
    for line in f:
        record = line[:6].strip()
        if record == "ATOM":
            chain   = line[21]
            try:
                res_num = int(line[22:26].strip())
            except ValueError:
                continue
            if chain == "A" and 66 <= res_num <= 404:
                receptor_lines.append(line)
        elif record == "HETATM":
            res_name = line[17:20].strip()
            if res_name == LIGAND_RESNAME:
                ligand_lines.append(line)
        elif record in ("HEADER", "TITLE", "REMARK"):
            receptor_lines.append(line)

receptor_lines.append("END\n")

with open(RECEPTOR_PDB, "w") as f:
    f.writelines(receptor_lines)

n_atoms = sum(1 for l in receptor_lines if l.startswith("ATOM"))
print(f"  {n_atoms} ATOM records written")
print(f"  {len(ligand_lines)} ligand ({LIGAND_RESNAME}) atoms found")

# ── 2. Compute binding site centroid ─────────────────────────────────────────

if not ligand_lines:
    print("ERROR: No ligand atoms found — cannot compute binding site")
    exit(1)

xs, ys, zs = [], [], []
for line in ligand_lines:
    try:
        xs.append(float(line[30:38]))
        ys.append(float(line[38:46]))
        zs.append(float(line[46:54]))
    except ValueError:
        continue

cx = sum(xs) / len(xs)
cy = sum(ys) / len(ys)
cz = sum(zs) / len(zs)

print(f"\nBinding site centroid ({LIGAND_RESNAME} centroid):")
print(f"  centre = ({cx:.2f}, {cy:.2f}, {cz:.2f})")
print(f"  box    = {BOX_SIZE} x {BOX_SIZE} x {BOX_SIZE} Angstrom")

with open("data/binding_site.txt", "w") as f:
    f.write(f"# Docking box parameters for 5-HT2A receptor (6WHA)\n")
    f.write(f"# Centre = centroid of co-crystallised ligand {LIGAND_RESNAME}\n")
    f.write(f"# Box size chosen to cover full orthosteric binding pocket\n")
    f.write(f"center_x = {cx:.3f}\n")
    f.write(f"center_y = {cy:.3f}\n")
    f.write(f"center_z = {cz:.3f}\n")
    f.write(f"size_x = {BOX_SIZE}\n")
    f.write(f"size_y = {BOX_SIZE}\n")
    f.write(f"size_z = {BOX_SIZE}\n")

print(f"  Saved data/binding_site.txt")

# ── 3. Convert to PDBQT ───────────────────────────────────────────────────────

print(f"\nConverting receptor to PDBQT (obabel, pH 7.4)...")

cmd = [
    "obabel",
    RECEPTOR_PDB,
    "-O", RECEPTOR_PDBQT,
    "-p", "7.4",
    "-xr",
]
result = subprocess.run(cmd, capture_output=True, text=True)

if result.stdout.strip():
    print(f"  {result.stdout.strip()}")
if result.stderr.strip():
    # obabel writes progress to stderr — show last line only
    last = [l for l in result.stderr.strip().split("\n") if l.strip()]
    if last:
        print(f"  {last[-1]}")

if os.path.exists(RECEPTOR_PDBQT):
    n_lines = sum(1 for _ in open(RECEPTOR_PDBQT))
    print(f"  Saved {RECEPTOR_PDBQT} ({n_lines} lines)")
else:
    print("  ERROR: PDBQT file not created")
    print("  Full stderr:")
    print(result.stderr)
    exit(1)

print("\nDone.")
print(f"  Receptor PDB:   {RECEPTOR_PDB}")
print(f"  Receptor PDBQT: {RECEPTOR_PDBQT}")
print(f"  Binding site:   data/binding_site.txt")
