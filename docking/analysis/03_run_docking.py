"""
Script 03: Run AutoDock Vina docking for all prepared Shulgin ligands.

Each ligand is docked in an isolated subprocess so that C++ internal_error
crashes in Vina don't abort the entire run. The main process collects
results from each subprocess and continues to the next ligand.

Outputs:
  results/docking_scores.csv
  results/poses/<compound>.pdbqt
  analysis/03_docking_report.txt

Usage:
  cd psychedelic-selectivity/docking/
  python analysis/03_run_docking.py
"""

import pandas as pd
import numpy as np
import os
import sys
import time
import subprocess
import json

os.makedirs("results/poses", exist_ok=True)

# ── Parameters ────────────────────────────────────────────────────────────────

RECEPTOR_PDBQT = "data/6WHA_receptor.pdbqt"
MANIFEST_CSV   = "data/ligands/ligand_manifest.csv"
BINDING_SITE   = "data/binding_site.txt"
EXHAUSTIVENESS = 8
N_POSES        = 9
ENERGY_RANGE   = 3.0
CPU            = 4

# ── Worker script (run per ligand in subprocess) ──────────────────────────────
# Written to a temp file and called as a subprocess.
# On success prints JSON to stdout. On failure exits with code 1.

WORKER_SCRIPT = """
import sys, json
from vina import Vina

receptor  = sys.argv[1]
ligand_in = sys.argv[2]
poses_out = sys.argv[3]
cx, cy, cz = float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6])
sx, sy, sz = float(sys.argv[7]), float(sys.argv[8]), float(sys.argv[9])
exhaustiveness = int(sys.argv[10])
n_poses        = int(sys.argv[11])
energy_range   = float(sys.argv[12])
cpu            = int(sys.argv[13])

v = Vina(sf_name="vina", cpu=cpu, verbosity=0)
v.set_receptor(receptor)
v.compute_vina_maps(center=[cx, cy, cz], box_size=[sx, sy, sz])
v.set_ligand_from_file(ligand_in)
v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
v.write_poses(poses_out, n_poses=n_poses, energy_range=energy_range, overwrite=True)
energies = v.energies(n_poses=n_poses, energy_range=energy_range)
result = {
    "best_score":   energies[0][0],
    "second_score": energies[1][0] if len(energies) > 1 else None,
    "n_poses":      len(energies),
}
print(json.dumps(result))
"""

WORKER_PATH = "/tmp/vina_worker.py"
with open(WORKER_PATH, "w") as f:
    f.write(WORKER_SCRIPT)

# ── Load binding site ─────────────────────────────────────────────────────────

box = {}
with open(BINDING_SITE) as f:
    for line in f:
        if line.startswith("#"): continue
        key, val = line.strip().split(" = ")
        box[key.strip()] = float(val.strip())

cx = box["center_x"]; cy = box["center_y"]; cz = box["center_z"]
sx = box["size_x"];   sy = box["size_y"];   sz = box["size_z"]

# ── Load manifest ─────────────────────────────────────────────────────────────

manifest = pd.read_csv(MANIFEST_CSV)
manifest = manifest[manifest["status"] == "OK"].reset_index(drop=True)
manifest = manifest.drop_duplicates(subset="pdbqt").reset_index(drop=True)
n_total  = len(manifest)

log = open("analysis/03_docking_report.txt", "w")
def p(*args, **kwargs):
    print(*args, **kwargs, file=log)

p("AutoDock Vina Docking — Shulgin compounds at 5-HT2A (6WHA)")
p("=" * 65)
p(f"Receptor:       {RECEPTOR_PDBQT}")
p(f"Box centre:     ({cx}, {cy}, {cz})")
p(f"Box size:       {sx} x {sy} x {sz} Angstrom")
p(f"Exhaustiveness: {EXHAUSTIVENESS}")
p(f"N poses:        {N_POSES}")
p(f"Energy range:   {ENERGY_RANGE} kcal/mol")
p(f"Ligands:        {n_total}")
p("")
p(f"{'#':>3} {'Compound':<28} {'Best score':>11} {'2nd pose':>9} {'Time(s)':>8} {'Status'}")
p("-" * 75)

# ── Dock each ligand in isolated subprocess ───────────────────────────────────

results = []

for i, row in manifest.iterrows():
    compound  = row["compound"]
    pdbqt_in  = row["pdbqt"]
    pdbqt_out = f"results/poses/{row['safe_name']}.pdbqt"

    if not os.path.exists(pdbqt_in):
        p(f"{i+1:>3} {compound:<28} {'n/a':>11} {'n/a':>9} {'n/a':>8} MISSING")
        results.append({"compound": compound, "best_score": np.nan,
                        "status": "missing"})
        continue

    t0 = time.time()

    cmd = [
        sys.executable, WORKER_PATH,
        RECEPTOR_PDBQT, pdbqt_in, pdbqt_out,
        str(cx), str(cy), str(cz),
        str(sx), str(sy), str(sz),
        str(EXHAUSTIVENESS), str(N_POSES), str(ENERGY_RANGE), str(CPU),
    ]

    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    elapsed = time.time() - t0

    if proc.returncode == 0 and proc.stdout.strip():
        try:
            res = json.loads(proc.stdout.strip())
            best   = res["best_score"]
            second = res["second_score"] if res["second_score"] else np.nan
            p(f"{i+1:>3} {compound:<28} {best:>11.3f} {second:>9.3f} "
              f"{elapsed:>8.1f} OK")
            print(f"  [{i+1}/{n_total}] {compound:<28} {best:.3f} kcal/mol  "
                  f"({elapsed:.0f}s)")
            results.append({"compound": compound, "best_score": best,
                            "second_score": second, "n_poses": res["n_poses"],
                            "elapsed_s": elapsed, "status": "OK",
                            "pdbqt_out": pdbqt_out})
        except Exception as e:
            p(f"{i+1:>3} {compound:<28} {'n/a':>11} {'n/a':>9} "
              f"{elapsed:>8.1f} JSON_ERROR: {e}")
            print(f"  [{i+1}/{n_total}] {compound:<28} JSON error: {e}")
            results.append({"compound": compound, "best_score": np.nan,
                            "status": f"json_error: {e}"})
    else:
        err = (proc.stderr or proc.stdout or "unknown")[:80].strip()
        p(f"{i+1:>3} {compound:<28} {'n/a':>11} {'n/a':>9} "
          f"{elapsed:>8.1f} FAILED: {err}")
        print(f"  [{i+1}/{n_total}] {compound:<28} FAILED: {err}")
        results.append({"compound": compound, "best_score": np.nan,
                        "status": f"failed: {err}"})

# ── Save and summarise ────────────────────────────────────────────────────────

results_df = pd.DataFrame(results)
results_df.to_csv("results/docking_scores.csv", index=False)

ok = results_df[results_df["status"] == "OK"]
p("")
p("=" * 65)
p("SUMMARY")
p("=" * 65)
p(f"Successfully docked: {len(ok)} / {len(results_df)}")
if len(ok) > 0:
    p(f"Score range: {ok['best_score'].min():.3f} to "
      f"{ok['best_score'].max():.3f} kcal/mol")
    p(f"\nTop 10 best binders:")
    for _, r in ok.nsmallest(10, "best_score").iterrows():
        p(f"  {r['compound']:<32} {r['best_score']:>8.3f} kcal/mol")
    p(f"\nBottom 5 weakest binders:")
    for _, r in ok.nlargest(5, "best_score").iterrows():
        p(f"  {r['compound']:<32} {r['best_score']:>8.3f} kcal/mol")

p(f"\nSaved results/docking_scores.csv")
log.close()

print(f"\nDone. {len(ok)}/{len(results_df)} docked successfully.")
print(f"Results: results/docking_scores.csv")
print(f"Report:  analysis/03_docking_report.txt")
