"""
AptaScope v2 — Branch 3 build step.

Reads public/data/aptascope_merged.json (13,306 records, already carries
predicted_mfe/mfe_structure/num_stems/num_loops from the merge step — no
folding needed here) and writes src/data/optimization_profiles.json: a
per-target-type statistical profile of what distinguishes top-25%-by-affinity
binders from the rest of the cohort, which the Analyzer's Optimizer panel
(src/utils/optimizer.js) uses at runtime to generate position-specific
suggestions for a user-pasted sequence.

Usage:
    python build_optimization_profiles.py
"""

import json
import math
from collections import Counter
from pathlib import Path

IN_PATH = Path("public/data/aptascope_merged.json")
OUT_PATH = Path("src/data/optimization_profiles.json")

# Relative-position bins (as a fraction of sequence length) rather than raw
# indices — records range from ~15nt to ~80nt, so "position 40" means
# something different in each. Binning by fraction lets a 20nt and an 80nt
# sequence both contribute to e.g. "3' end" preference.
NUM_BINS = 5
BIN_LABELS = ["5' end", "early", "middle", "late", "3' end"]
DINUCLEOTIDES = [a + b for a in "ATGC" for b in "ATGC"]
MIN_COHORT = 20  # below this, a target type falls back to the "All" profile


def percentile(sorted_vals, p):
    if not sorted_vals:
        return None
    if len(sorted_vals) == 1:
        return sorted_vals[0]
    k = (len(sorted_vals) - 1) * p
    f, c = math.floor(k), math.ceil(k)
    if f == c:
        return sorted_vals[int(k)]
    return sorted_vals[f] + (sorted_vals[c] - sorted_vals[f]) * (k - f)


def median(vals):
    s = sorted(vals)
    n = len(s)
    if n == 0:
        return None
    mid = n // 2
    return s[mid] if n % 2 else (s[mid - 1] + s[mid]) / 2


def norm_base(ch):
    return "T" if ch == "U" else ch


def position_bin(i, length):
    frac = i / length
    idx = min(NUM_BINS - 1, int(frac * NUM_BINS))
    return idx


def base_frequencies_by_bin(records):
    """{bin_idx: {base: freq}} across all positions of all records, pooled."""
    bin_counts = [Counter() for _ in range(NUM_BINS)]
    bin_totals = [0] * NUM_BINS
    for r in records:
        seq = r["sequence"]
        length = len(seq)
        for i, ch in enumerate(seq):
            b = position_bin(i, length)
            bin_counts[b][norm_base(ch)] += 1
            bin_totals[b] += 1
    out = []
    for b in range(NUM_BINS):
        total = bin_totals[b] or 1
        out.append({base: bin_counts[b].get(base, 0) / total for base in "ATGC"})
    return out


def dinucleotide_frequencies(records):
    counts = Counter()
    total = 0
    for r in records:
        seq = "".join(norm_base(c) for c in r["sequence"])
        for i in range(len(seq) - 1):
            counts[seq[i : i + 2]] += 1
            total += 1
    total = total or 1
    return {d: counts.get(d, 0) / total for d in DINUCLEOTIDES}


def build_profile(records, label):
    with_kd = [r for r in records if r["kd_nM"] is not None]
    with_kd.sort(key=lambda r: r["kd_nM"])
    cutoff = max(1, len(with_kd) // 4)
    top = with_kd[:cutoff]
    rest = with_kd[cutoff:] or with_kd

    gcs_top = sorted(r["gc_content"] for r in top)
    lengths_top = sorted(r["length"] for r in top)

    top_bins = base_frequencies_by_bin(top)
    rest_bins = base_frequencies_by_bin(rest)
    position_preferences = []
    for b in range(NUM_BINS):
        enrichment = {}
        for base in "ATGC":
            rest_f = rest_bins[b][base]
            top_f = top_bins[b][base]
            enrichment[base] = round(top_f / rest_f, 3) if rest_f > 0 else None
        position_preferences.append(
            {"label": BIN_LABELS[b], "top_freq": {k: round(v, 4) for k, v in top_bins[b].items()}, "enrichment": enrichment}
        )

    top_dinuc = dinucleotide_frequencies(top)
    rest_dinuc = dinucleotide_frequencies(rest)
    dinucleotide_enrichment = {
        d: round(math.log2(top_dinuc[d] / rest_dinuc[d]), 3) if top_dinuc[d] > 0 and rest_dinuc[d] > 0 else 0.0
        for d in DINUCLEOTIDES
    }

    top_with_struct = [r for r in top if r.get("predicted_mfe") is not None]
    rest_with_struct = [r for r in rest if r.get("predicted_mfe") is not None]
    mfe_per_nt_top = median([r["predicted_mfe"] / r["length"] for r in top_with_struct]) if top_with_struct else None
    mfe_per_nt_rest = median([r["predicted_mfe"] / r["length"] for r in rest_with_struct]) if rest_with_struct else None
    stems_per_nt_top = median([r["num_stems"] / r["length"] for r in top_with_struct]) if top_with_struct else None

    return {
        "label": label,
        "sample_size": len(with_kd),
        "top_cohort_size": len(top),
        "gc_range": [round(percentile(gcs_top, 0.25), 4), round(percentile(gcs_top, 0.75), 4)] if gcs_top else None,
        "length_range": [round(percentile(lengths_top, 0.25)), round(percentile(lengths_top, 0.75))] if lengths_top else None,
        "position_preferences": position_preferences,
        "dinucleotide_enrichment": dinucleotide_enrichment,
        "structural": {
            "optimal_mfe_per_nt": round(mfe_per_nt_top, 4) if mfe_per_nt_top is not None else None,
            "typical_mfe_per_nt": round(mfe_per_nt_rest, 4) if mfe_per_nt_rest is not None else None,
            "optimal_stems_per_nt": round(stems_per_nt_top, 4) if stems_per_nt_top is not None else None,
        },
    }


def main():
    with open(IN_PATH) as f:
        records = json.load(f)

    profiles = {"All": build_profile(records, "All targets")}

    by_type = {}
    for r in records:
        by_type.setdefault(r["target_type"], []).append(r)

    for target_type, recs in by_type.items():
        cohort_with_kd = sum(1 for r in recs if r["kd_nM"] is not None)
        if cohort_with_kd < MIN_COHORT:
            continue  # too few affinity-labeled records to fit a stable profile
        profiles[target_type] = build_profile(recs, target_type)

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_PATH, "w") as f:
        json.dump(profiles, f, indent=2)

    print(f"Wrote {len(profiles)} profiles to {OUT_PATH}:")
    for name, p in profiles.items():
        print(f"  {name}: n={p['sample_size']}, gc_range={p['gc_range']}, length_range={p['length_range']}")


if __name__ == "__main__":
    main()
