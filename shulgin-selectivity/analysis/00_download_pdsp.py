import requests
import csv
import time
import sys

BASE_URL = "https://pdspdb.unc.edu/databases/kiDownload/search.php"
OUT_FILE = "analysis/KiDatabase.csv"
PAGE_SIZE = 20
MAX_RECORDS = 100000

FIELDS = ["number", "name", "unigene", "ligandid", "ligandname",
          "smiles", "cas", "nsc", "hotligand", "species",
          "source", "kinote", "kival", "reference", "link"]

print("Downloading PDSP Ki database...")
total = 0
offset = 0
errors = 0

with open(OUT_FILE, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=FIELDS, extrasaction="ignore")
    writer.writeheader()
    while offset < MAX_RECORDS:
        try:
            r = requests.get(BASE_URL, params={"offset": offset}, timeout=30)
            r.raise_for_status()
            records = r.json().get("records", [])
            if not records:
                print(f"\nNo more records at offset {offset}. Done.")
                break
            writer.writerows(records)
            f.flush()
            total += len(records)
            if offset % 1000 == 0:
                print(f"  offset={offset:>6}  total={total:>6}", end="\r")
                sys.stdout.flush()
            offset += PAGE_SIZE
            time.sleep(0.05)
        except Exception as e:
            errors += 1
            print(f"\n  Error at offset {offset}: {e}")
            if errors > 10:
                print("Too many errors, stopping.")
                break
            time.sleep(2)

print(f"\nDone. {total} records, {errors} errors.")

import pandas as pd
df = pd.read_csv(OUT_FILE, low_memory=False)
print(f"Shape: {df.shape}")
print(df["name"].value_counts().head(10).to_string())
