"""
run_all_seeds.py

Runs all iSIC datasets multiple times with different random seeds,
then produces per-dataset statistics tables: min/avg/max of the clustering
function value, average CPU time, average D-B index, and average Dunn index,
for every number of clusters nc.

Usage:
    cd src/
    python3 run_all_seeds.py [nroseeds]

    nroseeds   (optional) number of seeds to use, default 10.

The seed sequence is fixed and reproducible: always the same seeds are used
in the same order regardless of how many times the script is run.
"""

import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Number of seeds to use.  Can also be overridden from the command line:
#   python3 run_all_seeds.py 5
NROSEEDS = int(sys.argv[1]) if len(sys.argv) > 1 else 10

# Fixed, reproducible seed sequence.  Always use the first NROSEEDS of these.
ALL_SEEDS = [
    12345, 54321, 11111, 22222, 33333,
    44444, 55555, 66666, 77777, 88888,
    13579, 24680, 98765, 86420, 11223,
    33445, 55667, 77889, 19283, 47561,
]
SEEDS = ALL_SEEDS[:NROSEEDS]

SRC_DIR = Path(__file__).parent.resolve()
INIT_FILE = SRC_DIR / "initisic.f95"
BACKUP_FILE = SRC_DIR / "initisic_backup_seeds.f95"
SUMMARY_OUT = SRC_DIR / "isic_all_seeds_summary.txt"

DATA_BASE = "../../../data/data_for_clustering/large_data_for_clustering"

DATASETS = [
    {"id": 11, "file": "isolet.txt",          "nft": 617,  "nrecord":   7797},
    {"id": 14, "file": "gisette.txt",         "nft": 5000, "nrecord":  13500},
    {"id":  1, "file": "drift.txt",           "nft": 128,  "nrecord":  13910},
    {"id":  2, "file": "EEG_Eye_State2.txt",  "nft":  14,  "nrecord":  14980},
    {"id":  3, "file": "d15112.txt",          "nft":   2,  "nrecord":  15112},
    {"id": 12, "file": "OnlineNews.txt",      "nft":  58,  "nrecord":  39644},
    {"id":  4, "file": "KEGG2.txt",           "nft":  20,  "nrecord":  53413},
    {"id":  5, "file": "shuttle.txt",         "nft":   9,  "nrecord":  58000},
    {"id":  9, "file": "Sensorless.txt",      "nft":  48,  "nrecord":  58509},
    {"id": 15, "file": "mfcc.txt",            "nft":  58,  "nrecord":  85134},
    {"id":  6, "file": "pla85900.txt",        "nft":   2,  "nrecord":  85900},
    {"id": 10, "file": "music.txt",         "nft":  518,  "nrecord": 106574},
    {"id": 13, "file": "MiniBooNE_PID.txt",   "nft":  50,  "nrecord": 130064},
    {"id": 16, "file": "prot_homo.txt",       "nft":  74,  "nrecord": 145751},
    {"id": 17, "file": "rq3.txt",             "nft":   7,  "nrecord": 200000},
    {"id":  7, "file": "Skin_NonSkin.txt",    "nft":   3,  "nrecord": 245057},
    {"id":  8, "file": "3D.txt",              "nft":   3,  "nrecord": 434874},
    {"id": 18, "file": "Covertype.txt",       "nft":  10,  "nrecord": 581012},
    {"id": 20, "file": "census.txt",       "nft":  68,  "nrecord": 2458285},
    {"id": 21, "file": "bitcoin.txt",       "nft":  8,  "nrecord": 2916697},
    {"id": 22, "file": "hepmass.txt",       "nft":  28,  "nrecord": 10500000},
    #{"id": 19, "file": "cord.csv",         "nft":  768,  "nrecord": 1056660},
]

# ---------------------------------------------------------------------------
# Regex patterns — same as run_all_datasets.py
# ---------------------------------------------------------------------------

RE_INFILE = re.compile(
    r"^(\s+infile\s*=\s*')([^']+)(',\s*&.*)",
    re.MULTILINE,
)
RE_NFT = re.compile(
    r"^(\s+nft\s*=\s*)(\d+[-\d]*)(,\s*&.*)",
    re.MULTILINE,
)
RE_NRECORD = re.compile(
    r"^(\s+nrecord\s*=\s*)(\d+)(,\s*&.*)",
    re.MULTILINE,
)
RE_NREC_USED = re.compile(
    r"^(\s+nrec_used\s*=\s*)([^,\n]+)(,\s*&.*)",
    re.MULTILINE,
)
RE_OUTFILE1 = re.compile(
    r"^(\s+outfile1\s*=\s*')([^']+)('.*)",
    re.MULTILINE,
)


def patch_init_file(ds: dict) -> str:
    """Patch initisic.f95 for this dataset.  Returns outfile1 name."""
    text = INIT_FILE.read_text()

    infile_path = f"{DATA_BASE}/{ds['file']}"
    outfile1_name = "isic_results_current.txt"

    text, n = RE_INFILE.subn(lambda m: m.group(1) + infile_path + m.group(3), text)
    if n != 1:
        raise ValueError(f"Expected 1 active infile line, found {n}")

    text, n = RE_NFT.subn(lambda m: m.group(1) + str(ds["nft"]) + m.group(3), text)
    if n != 1:
        raise ValueError(f"Expected 1 active nft line, found {n}")

    text, n = RE_NRECORD.subn(lambda m: m.group(1) + str(ds["nrecord"]) + m.group(3), text)
    if n != 1:
        raise ValueError(f"Expected 1 active nrecord line, found {n}")

    text, n = RE_NREC_USED.subn(lambda m: m.group(1) + "nrecord" + m.group(3), text)
    if n != 1:
        raise ValueError(f"Expected 1 active nrec_used line, found {n}")

    text, n = RE_OUTFILE1.subn(lambda m: m.group(1) + outfile1_name + m.group(3), text)
    if n != 1:
        raise ValueError(f"Expected 1 active outfile1 line, found {n}")

    INIT_FILE.write_text(text)
    return outfile1_name


def parse_results(path: Path) -> dict:
    """Parse an outfile1 result file.

    Returns {nc: {"f": float, "db": float, "dn": float, "cpu": float}}
    for every nc row found.  Non-data lines (headers, blanks) are skipped.
    """
    results = {}
    for line in path.read_text().splitlines():
        tokens = line.split()
        if len(tokens) < 6:
            continue
        try:
            nc = int(tokens[0])
            f = float(tokens[1])
            db = float(tokens[2])
            dn = float(tokens[3])
            # tokens[4] = tnorm, tokens[5] = cpu_time
            cpu = float(tokens[5])
        except ValueError:
            continue
        results[nc] = {"f": f, "db": db, "dn": dn, "cpu": cpu}
    return results


def aggregate(all_runs: list[dict]) -> dict:
    """Given a list of per-seed result dicts, compute statistics per nc.

    Returns {nc: {"f_min", "f_avg", "f_max", "cpu_avg", "db_avg", "dn_avg"}}
    Only nc values that appear in ALL runs are included.
    """
    # Collect all nc values present in every run
    common_ncs = set(all_runs[0].keys())
    for run in all_runs[1:]:
        common_ncs &= set(run.keys())

    stats = {}
    for nc in sorted(common_ncs):
        f_vals   = [r[nc]["f"]   for r in all_runs]
        db_vals  = [r[nc]["db"]  for r in all_runs]
        dn_vals  = [r[nc]["dn"]  for r in all_runs]
        cpu_vals = [r[nc]["cpu"] for r in all_runs]
        n = len(all_runs)
        stats[nc] = {
            "f_min":   min(f_vals),
            "f_avg":   sum(f_vals)   / n,
            "f_max":   max(f_vals),
            "cpu_avg": sum(cpu_vals) / n,
            "db_avg":  sum(db_vals)  / n,
            "dn_avg":  sum(dn_vals)  / n,
            "n_runs":  n,
        }
    return stats


def format_stats_table(stats: dict) -> str:
    """Format the statistics dict as a plain-text table."""
    header = (
        f"{'nc':>4} | {'f_min':>18} | {'f_avg':>18} | {'f_max':>18} "
        f"| {'cpu_avg':>9} | {'db_avg':>9} | {'dn_avg':>9} | runs"
    )
    sep = "-" * len(header)
    lines = [header, sep]
    for nc, s in sorted(stats.items()):
        lines.append(
            f"{nc:>4} | {s['f_min']:>18.4f} | {s['f_avg']:>18.4f} | {s['f_max']:>18.4f} "
            f"| {s['cpu_avg']:>9.2f} | {s['db_avg']:>9.4f} | {s['dn_avg']:>9.4f} | {s['n_runs']}"
        )
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print(f"iSIC multi-seed test runner")
    print(f"  Datasets : {len(DATASETS)}")
    print(f"  Seeds    : {NROSEEDS}  {SEEDS}")
    print(f"  Output   : {SUMMARY_OUT.name}")

    print(f"\nBacking up {INIT_FILE.name} -> {BACKUP_FILE.name}")
    shutil.copy2(INIT_FILE, BACKUP_FILE)

    with SUMMARY_OUT.open("w") as out:
        out.write("iSIC — multi-seed statistics summary\n")
        out.write(f"Seeds used ({NROSEEDS}): {SEEDS}\n")
        out.write("=" * 80 + "\n\n")

        for ds in DATASETS:
            ds_tag = f"Dataset {ds['id']}: {ds['file']}  (nft={ds['nft']}, nrecord={ds['nrecord']})"
            print(f"\n{'='*60}\n  {ds_tag}\n{'='*60}")

            # --- Patch and compile once for this dataset ---
            try:
                patch_init_file(ds)
            except ValueError as e:
                msg = f"[SKIP] Regex patch failed: {e}"
                print(msg)
                out.write(f"\n{ds_tag}\n{msg}\n")
                shutil.copy2(BACKUP_FILE, INIT_FILE)
                continue

            print("  Compiling...")
            result = subprocess.run(["make"], cwd=SRC_DIR, capture_output=True, text=True)
            if result.returncode != 0:
                msg = f"[SKIP] make failed:\n{result.stderr}"
                print(msg)
                out.write(f"\n{ds_tag}\n{msg}\n")
                shutil.copy2(BACKUP_FILE, INIT_FILE)
                continue
            print("  Compilation successful.")

            # --- Run once per seed ---
            per_ds_file = SRC_DIR / "isic_results_current.txt"
            all_runs = []
            for seed in SEEDS:
                print(f"  Seed {seed} ...", end="", flush=True)

                run_result = subprocess.run(
                    ["./isic", str(seed)],
                    cwd=SRC_DIR,
                    capture_output=True,
                    text=True,
                )
                if run_result.returncode != 0:
                    print(f" [FAIL] {run_result.stderr[:120]}")
                    continue
                if not per_ds_file.exists():
                    print(" [FAIL] output file missing")
                    continue

                run_data = parse_results(per_ds_file)
                if not run_data:
                    print(" [FAIL] no data rows parsed")
                    continue

                all_runs.append(run_data)
                print(f" done ({len(run_data)} nc values)")

            # Restore backup before next dataset
            shutil.copy2(BACKUP_FILE, INIT_FILE)

            # --- Write statistics ---
            out.write(f"\n{'='*80}\n")
            out.write(f"  {ds_tag}\n")
            out.write(f"  Seeds: {SEEDS[:len(all_runs)]}  ({len(all_runs)} successful runs)\n")
            out.write(f"{'='*80}\n")

            if not all_runs:
                out.write("  [No successful runs.]\n")
                continue

            stats = aggregate(all_runs)
            out.write(format_stats_table(stats))
            out.write("\n")
            out.flush()

            print(f"  Statistics written for {len(stats)} nc values.")

    print(f"\nAll done. Summary: {SUMMARY_OUT}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nInterrupted. Restoring initisic.f95 ...")
        shutil.copy2(BACKUP_FILE, INIT_FILE)
        sys.exit(1)
    finally:
        if BACKUP_FILE.exists():
            shutil.copy2(BACKUP_FILE, INIT_FILE)
