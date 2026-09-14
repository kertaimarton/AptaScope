"""
AptaScope — DNA structure bug fix (one-off re-fold).

merge_databases.py's structure-filling step converted every sequence's
T -> U and ran ViennaRNA's default RNA fold on it, including DNA aptamers.
The RNA model allows G-U "wobble" pairs, which come back out as G-T once
remapped to DNA letters — a real tolerance in RNA, not a real DNA base
pair. Across the merged dataset, 68-81% of DNA records (depending on
source) carried at least one such non-canonical pair, including records
whose structure/MFE came straight from AptaNexus's own source data (not
just the ones this pipeline folded itself) — see predict_structure_if_missing
in merge_databases.py, now fixed for future re-runs.

This script re-folds every DNA record already in
public/data/aptascope_merged.json using ViennaRNA's actual DNA parameter set
(Mathews 2004) with wobble pairing disabled (noGU=1), overwrites
predicted_mfe/mfe_structure/num_stems/num_loops, and regenerates
dataset_stats.json. RNA records are untouched. Run this once; going
forward, merge_databases.py does this correctly itself.

Usage:
    python refold_dna_canonical.py
"""

import json
import re
from pathlib import Path

import RNA

from merge_databases import compute_stats

MERGED_PATH = Path("public/data/aptascope_merged.json")
STATS_PATH = Path("src/data/dataset_stats.json")


def count_stems_loops(dot_bracket):
    if not dot_bracket:
        return None, None
    stems = dot_bracket.count("(")
    loops = len(re.findall(r"\(\.+\)", dot_bracket))
    return stems, loops


def main():
    with open(MERGED_PATH) as f:
        records = json.load(f)

    md = RNA.md()
    md.noGU = 1
    RNA.params_load_DNA_Mathews2004()

    dna_records = [r for r in records if r["aptamer_type"] == "DNA"]
    print(f"Re-folding {len(dna_records)} DNA records with DNA-specific parameters (noGU)...")

    changed = 0
    skipped_short = 0
    for r in dna_records:
        seq = r["sequence"]
        if len(seq) < 3:
            r["predicted_mfe"], r["mfe_structure"], r["num_stems"], r["num_loops"] = None, None, None, None
            skipped_short += 1
            continue
        fc = RNA.fold_compound(seq, md)
        structure, mfe = fc.mfe()
        if structure != r.get("mfe_structure"):
            changed += 1
        r["predicted_mfe"] = round(float(mfe), 2)
        r["mfe_structure"] = structure
        r["num_stems"], r["num_loops"] = count_stems_loops(structure)

    print(f"{changed}/{len(dna_records)} DNA records had their structure changed ({skipped_short} too short to fold).")

    with open(MERGED_PATH, "w") as f:
        json.dump(records, f, separators=(",", ":"))
    print(f"Wrote {len(records)} records to {MERGED_PATH}")

    stats = compute_stats(records)
    with open(STATS_PATH, "w") as f:
        json.dump(stats, f, indent=2)
    print(f"Wrote stats to {STATS_PATH}")
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
