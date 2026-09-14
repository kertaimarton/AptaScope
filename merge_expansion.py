"""
AptaScope — merge the aptascope_datamore/ data expansion into the live
merged dataset.

Folds two pools into public/data/aptascope_merged.json:
  - src/data/aptamer_data.json (1,757 UTexas/Ribocentre records, pulled
    earlier and never wired in)
  - aptascope_datamore/data_pulls/aptascope_expansion.json (5,758 records
    from AptBacterialDB, UTexas, AptCancerDB, Aptagen, Ribocentre, RNAapt3D)

Deduped by SEQUENCE alone (not sequence+target_name) against the live
13,306-record dataset and against each other: measured overlap shows the
vast majority of collisions are the same physical aptamer reported under a
different target-name spelling across databases (e.g. "CEA" vs
"Carcinoembryonic antigen (CEA)"), not distinct aptamers that happen to
share a sequence — a (sequence, target_name) key would fail to catch most
of these and bloat the dataset with near-duplicate rows for one molecule.

Reuses merge_databases.py's cleaning/classification/folding/feature
functions rather than reimplementing them, so DNA structures come out with
the same corrected DNA-only base-pairing rules (no G-T wobble) already
applied to the rest of the dataset this session.
"""

import json
import re
from pathlib import Path

from merge_databases import (
    clean_sequence,
    detect_aptamer_type,
    clean_target_name,
    canonicalize_target,
    classify_aptanexus_target_type,
    predict_structure_if_missing,
    compute_features,
    compute_stats,
    completeness_score,
    _plausible_kd,
)

EXPANSION_PATH = Path("aptascope_datamore/data_pulls/aptascope_expansion.json")
LEGACY_PATH = Path("src/data/aptamer_data.json")
MERGED_PATH = Path("public/data/aptascope_merged.json")
STATS_PATH = Path("src/data/dataset_stats.json")

# Two sources' raw target_type labels use words classify_aptanexus_target_type
# doesn't recognize; remap to its vocabulary before classifying.
RAW_TARGET_TYPE_ALIASES = {
    "small organic": "small molecule",
    "antibody": "protein",
}

# AptCancerDB/AptBacterialDB records targeting a whole cell/organism (not a
# single protein) despite both sources hardcoding target_type="protein" for
# every record regardless of content (confirmed in pull_all_sources.py's
# normalize_aptcancerdb/normalize_aptbacterialdb — a constant, not derived).
WHOLE_ORGANISM_RE = re.compile(r"whole[\s-]?(cell|organism|bacteri)", re.I)

# Well-known targets (canonicalize_target's own alias table) whose real type
# is unambiguous — lets utexas/rnaapt3d records with no target_type field at
# all still classify correctly when the name resolves to one of these.
KNOWN_ALIAS_TYPES = {
    "Thrombin": "Protein", "VEGF": "Protein", "PDGF-AA": "Protein",
    "PDGF-BB": "Protein", "PDGF-AB": "Protein",
    "HIV-1 Reverse Transcriptase": "Protein", "HIV-1 Rev": "Protein",
    "Taq DNA Polymerase": "Protein", "IgE": "Protein", "bFGF": "Protein",
    "ATP": "Small Molecule",
}

CANONICAL_DEFAULTS = {
    "pkd": None, "quality_tier": None, "journal": None, "article_title": None,
    "external_id": None, "external_name": None, "gene_symbol": None,
    "id_type": None, "buffer": None, "is_best": None, "sequence_id": None,
    "protein_name": None, "pubmed_id": None,
    "cancer_type": None, "bacterial_species": None,
}

DOI_PATTERN = re.compile(r"^10\.")


def classify_expansion_target_type(name, raw_target_type, source, uniprot_id):
    if WHOLE_ORGANISM_RE.search(name or ""):
        return "Microorganism" if source == "aptbacterialdb" else "Cell"

    if source not in ("aptcancerdb", "aptbacterialdb"):
        # Only these two sources hardcode a meaningless constant; everyone
        # else's raw label is real per-record signal — trust it first.
        raw = (raw_target_type or "").strip().lower()
        raw = RAW_TARGET_TYPE_ALIASES.get(raw, raw)
        mapped = classify_aptanexus_target_type(raw)
        if mapped != "Other":
            return mapped

    if uniprot_id:
        return "Protein"
    if name in KNOWN_ALIAS_TYPES:
        return KNOWN_ALIAS_TYPES[name]
    return classify_aptanexus_target_type(name or "")


def split_doi_or_pubmed(raw_doi):
    if not isinstance(raw_doi, str):
        return None, None
    s = raw_doi.strip()
    if not s:
        return None, None
    if DOI_PATTERN.match(s):
        return s, None
    if s.isdigit():
        # Several sources mislabel a bare PubMed ID as "doi" — don't build a
        # doi.org link out of it.
        return None, s
    return None, None


def normalize_record(raw):
    sequence = clean_sequence(raw.get("sequence"))
    if sequence is None:
        return None

    aptamer_type = detect_aptamer_type(sequence, raw.get("type") or raw.get("aptamer_type") or "")
    target_name = canonicalize_target(clean_target_name(raw.get("target_name")))
    source = raw.get("source")
    uniprot_id = raw.get("uniprot_id") or None
    target_type = classify_expansion_target_type(target_name, raw.get("target_type"), source, uniprot_id)

    # Only trust a genuine supplied fold (rnaapt3d's dot_bracket). Everything
    # else — including a bare leftover predicted_mfe number with no
    # structure string behind it (src/data/aptamer_data.json) — gets a fresh
    # fold from predict_structure_if_missing() below, so an MFE value is
    # never detached from the structure it's supposed to describe.
    mfe_structure = raw.get("dot_bracket") or None

    kd_nM = raw.get("kd_nM")
    kd_nM = _plausible_kd(float(kd_nM)) if isinstance(kd_nM, (int, float)) else None

    year = raw.get("year")
    year = int(year) if isinstance(year, (int, float)) else None

    doi, pubmed_id = split_doi_or_pubmed(raw.get("doi"))

    record = {
        "sequence": sequence,
        "aptamer_type": aptamer_type,
        "length": len(sequence),
        "target_name": target_name,
        "target_type": target_type,
        "kd_nM": kd_nM,
        "predicted_mfe": None,
        "mfe_structure": mfe_structure,
        "year": year,
        "doi": doi,
        "pubmed_id": pubmed_id,
        "selex_method": raw.get("selex_method") or None,
        "uniprot_id": uniprot_id,
        "sequence_id": raw.get("aptamer_name") or None,
        "cancer_type": raw.get("cancer_type") or None,
        "bacterial_species": raw.get("bacterial_species") or None,
        "source": source,
    }
    for key, default in CANONICAL_DEFAULTS.items():
        record.setdefault(key, default)
    return record


def load_pool(path):
    with open(path) as f:
        rows = json.load(f)
    out = []
    dropped_short = 0
    dropped_invalid = 0
    for row in rows:
        if row.get("quality_flag") == "short_sequence":
            dropped_short += 1
            continue
        rec = normalize_record(row)
        if rec is None:
            dropped_invalid += 1
            continue
        out.append(rec)
    print(f"[{path.name}] loaded {len(out)}, dropped short_sequence={dropped_short} invalid={dropped_invalid}")
    return out


def dedupe_by_sequence(new_records, existing_sequences):
    best = {}
    for r in new_records:
        seq = r["sequence"]
        if seq in existing_sequences:
            continue
        cur = best.get(seq)
        if cur is None:
            best[seq] = r
        elif completeness_score(r) > completeness_score(cur):
            if cur["target_name"] != r["target_name"]:
                print(f"[dedupe] {seq[:24]}...: keeping {r['target_name']!r} over {cur['target_name']!r}")
            best[seq] = r
    return list(best.values())


def enrich_and_id(records):
    counters = {}
    enriched = []
    for r in records:
        predicted_mfe, mfe_structure = predict_structure_if_missing(
            r["sequence"], r["aptamer_type"], r["predicted_mfe"], r["mfe_structure"]
        )
        r["predicted_mfe"] = predicted_mfe
        r["mfe_structure"] = mfe_structure
        features = compute_features(r["sequence"], r["aptamer_type"], mfe_structure)

        n = counters.get(r["source"], 0) + 1
        counters[r["source"]] = n
        full_record = {"id": f"{r['source']}_{n:05d}", **r, **features}
        enriched.append(full_record)
        if len(enriched) % 500 == 0:
            print(f"  enriched {len(enriched)}/{len(records)}")
    return enriched


def main():
    with open(MERGED_PATH) as f:
        current_merged = json.load(f)
    existing_sequences = {r["sequence"] for r in current_merged}
    print(f"Existing merged dataset: {len(current_merged)} records")

    legacy = load_pool(LEGACY_PATH)
    expansion = load_pool(EXPANSION_PATH)
    combined = legacy + expansion
    print(f"Combined candidate pool (legacy+expansion): {len(combined)}")

    net_new = dedupe_by_sequence(combined, existing_sequences)
    print(f"Net-new records after sequence dedup: {len(net_new)}")

    enriched = enrich_and_id(net_new)

    all_records = current_merged + enriched
    with open(MERGED_PATH, "w") as f:
        json.dump(all_records, f, separators=(",", ":"))
    print(f"Wrote {len(all_records)} total records to {MERGED_PATH}")

    stats = compute_stats(all_records)
    with open(STATS_PATH, "w") as f:
        json.dump(stats, f, indent=2)
    print(f"Wrote stats to {STATS_PATH}")
    print(json.dumps(stats, indent=2))

    if LEGACY_PATH.exists():
        LEGACY_PATH.unlink()
        print(f"Deleted {LEGACY_PATH} (superseded by this merge)")


if __name__ == "__main__":
    main()
