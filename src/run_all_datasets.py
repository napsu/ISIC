"""
run_all_datasets.py

Loops over all ISIC datasets, patches initisic.f95 for each,
recompiles, runs isic, and appends each summary to one combined file.

Usage:
    cd src/
    python3 run_all_datasets.py
"""

import re
import shutil
import subprocess
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

SRC_DIR = Path(__file__).parent.resolve()
INIT_FILE = SRC_DIR / "initisic.f95"
BACKUP_FILE = SRC_DIR / "initisic_backup.f95"
COMBINED_OUT = SRC_DIR / "isic_all_results.txt"

DATA_BASE = "../../../data/data_for_clustering/large_data_for_clustering"

DATASETS = [
    {"id":  1, "file": "drift.txt",          "nft": 128,  "nrecord":  13910},
    {"id":  2, "file": "EEG_Eye_State2.txt",  "nft":  14,  "nrecord":  14980},
    {"id":  3, "file": "d15112.txt",          "nft":   2,  "nrecord":  15112},
    {"id":  4, "file": "KEGG2.txt",           "nft":  20,  "nrecord":  53413},
    {"id":  5, "file": "shuttle.txt",         "nft":   9,  "nrecord":  58000},
    {"id":  6, "file": "pla85900.txt",        "nft":   2,  "nrecord":  85900},
    {"id":  7, "file": "Skin_NonSkin.txt",    "nft":   3,  "nrecord": 245057},
    {"id":  8, "file": "3D.txt",              "nft":   3,  "nrecord": 434874},
    {"id":  9, "file": "Sensorless.txt",      "nft":  48,  "nrecord":  58509},
    {"id": 11, "file": "isolet.txt",          "nft": 617,  "nrecord":   7797},
    {"id": 12, "file": "OnlineNews.txt",      "nft":  58,  "nrecord":  39644},
    {"id": 13, "file": "MiniBooNE_PID.txt",   "nft":  50,  "nrecord": 130064},
    {"id": 14, "file": "gisette.txt",         "nft": 5000, "nrecord":  13500},
]

# ---------------------------------------------------------------------------
# Regex patterns for the four active (non-commented) lines to replace
# ---------------------------------------------------------------------------

# Matches the uncommented infile line, e.g.:
#     infile = '../../../.../foo.txt',&      ! ...
RE_INFILE = re.compile(
    r"^(\s+infile\s*=\s*')([^']+)(',\s*&.*)",
    re.MULTILINE,
)

# Matches the uncommented nft line, e.g.:
#     nft = 5000, &               !14
RE_NFT = re.compile(
    r"^(\s+nft\s*=\s*)(\d+[-\d]*)(,\s*&.*)",
    re.MULTILINE,
)

# Matches the uncommented nrecord line, e.g.:
#     nrecord = 13500, &      !14
RE_NRECORD = re.compile(
    r"^(\s+nrecord\s*=\s*)(\d+)(,\s*&.*)",
    re.MULTILINE,
)

# Matches the uncommented nrec_used line — either form:
#     nrec_used = nrecord, &  ...
#     nrec_used = int(MAX(...)), &  ...
RE_NREC_USED = re.compile(
    r"^(\s+nrec_used\s*=\s*)([^,\n]+)(,\s*&.*)",
    re.MULTILINE,
)

# Matches the uncommented outfile1 line, e.g.:
#     outfile1 = 'bigclust_ncand1_indices_gisette.txt'      ! ...
RE_OUTFILE1 = re.compile(
    r"^(\s+outfile1\s*=\s*')([^']+)('.*)",
    re.MULTILINE,
)


def patch_init_file(ds: dict) -> str:
    """Return the name of the per-dataset outfile1 (written by isic)."""
    text = INIT_FILE.read_text()

    infile_path = f"{DATA_BASE}/{ds['file']}"
    outfile1_name = f"isic_results_{ds['file'].replace('.txt', '')}.txt"

    def replace_infile(m):
        return m.group(1) + infile_path + m.group(3)

    def replace_nft(m):
        return m.group(1) + str(ds["nft"]) + m.group(3)

    def replace_nrecord(m):
        return m.group(1) + str(ds["nrecord"]) + m.group(3)

    def replace_nrec_used(m):
        return m.group(1) + "nrecord" + m.group(3)

    def replace_outfile1(m):
        return m.group(1) + outfile1_name + m.group(3)

    text, n = RE_INFILE.subn(replace_infile, text)
    if n != 1:
        raise ValueError(f"Expected 1 active infile line, found {n}")

    text, n = RE_NFT.subn(replace_nft, text)
    if n != 1:
        raise ValueError(f"Expected 1 active nft line, found {n}")

    text, n = RE_NRECORD.subn(replace_nrecord, text)
    if n != 1:
        raise ValueError(f"Expected 1 active nrecord line, found {n}")

    text, n = RE_NREC_USED.subn(replace_nrec_used, text)
    if n != 1:
        raise ValueError(f"Expected 1 active nrec_used line, found {n}")

    text, n = RE_OUTFILE1.subn(replace_outfile1, text)
    if n != 1:
        raise ValueError(f"Expected 1 active outfile1 line, found {n}")

    INIT_FILE.write_text(text)
    return outfile1_name


# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------

def main():
    print(f"Backing up {INIT_FILE.name} -> {BACKUP_FILE.name}")
    shutil.copy2(INIT_FILE, BACKUP_FILE)

    with COMBINED_OUT.open("w") as combined:
        combined.write("ISIC — results for all datasets\n")
        combined.write("=" * 70 + "\n\n")

        for ds in DATASETS:
            tag = f"Dataset {ds['id']}: {ds['file']}  (nft={ds['nft']}, nrecord={ds['nrecord']})"
            print(f"\n{'='*60}")
            print(f"  {tag}")
            print(f"{'='*60}")

            try:
                outfile1_name = patch_init_file(ds)
            except ValueError as e:
                msg = f"[SKIP] Regex patch failed for {ds['file']}: {e}\n"
                print(msg)
                combined.write(f"\n--- {tag} ---\n{msg}\n")
                shutil.copy2(BACKUP_FILE, INIT_FILE)
                continue

            # Compile
            print("  Compiling...")
            result = subprocess.run(
                ["make"],
                cwd=SRC_DIR,
                capture_output=True,
                text=True,
            )
            if result.returncode != 0:
                msg = f"[SKIP] make failed:\n{result.stderr}\n"
                print(msg)
                combined.write(f"\n--- {tag} ---\n{msg}\n")
                shutil.copy2(BACKUP_FILE, INIT_FILE)
                continue
            print("  Compilation successful.")

            # Run
            print("  Running isic...")
            result = subprocess.run(
                ["./isic"],
                cwd=SRC_DIR,
                capture_output=True,
                text=True,
            )
            if result.returncode != 0:
                msg = f"[SKIP] isic failed:\n{result.stderr}\n{result.stdout}\n"
                print(msg)
                combined.write(f"\n--- {tag} ---\n{msg}\n")
                shutil.copy2(BACKUP_FILE, INIT_FILE)
                continue
            print("  Run complete.")

            # Collect per-dataset summary
            per_ds_file = SRC_DIR / outfile1_name
            combined.write(f"\n{'='*70}\n")
            combined.write(f"  {tag}\n")
            combined.write(f"{'='*70}\n")
            if per_ds_file.exists():
                combined.write(per_ds_file.read_text())
            else:
                combined.write(f"[WARNING] Output file {outfile1_name} not found.\n")
            combined.flush()

            # Restore backup so the next patch starts from a clean state
            shutil.copy2(BACKUP_FILE, INIT_FILE)

    print(f"\nAll datasets processed. Combined results: {COMBINED_OUT}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nInterrupted. Restoring original initisic.f95 ...")
        shutil.copy2(BACKUP_FILE, INIT_FILE)
        sys.exit(1)
    finally:
        # Always restore the backup, even on unexpected errors
        if BACKUP_FILE.exists():
            shutil.copy2(BACKUP_FILE, INIT_FILE)
