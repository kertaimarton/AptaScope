"""
AptaScope v2 — Phase 0 database merge.

Merges AptaNexus (aptanexus_raw.jsonl, ~12,535 records) and AptaDB
(aptadb_raw.json, ~1,350 records) into a single unified schema per
APTASCOPE_V2_SPEC.md, and writes src/data/aptascope_merged.json.

Usage:
    python merge_databases.py \\
        --aptanexus aptascope_aptanexus/aptanexus_raw.jsonl \\
        --aptadb aptascope_aptadb/aptadb_raw.json \\
        --output src/data/aptascope_merged.json

Notes on data quality (discovered while building this):
  - AptaDB's raw `sequence` field is already clean ATGCU-only text for all
    1,350 records. Its own precomputed `length`/`gc_content` fields are
    NOT trustworthy though (e.g. Apta_1's `length` reports 106 while its
    actual sequence is 37 nt) — so gc_content/frequencies/complexity/
    g-quadruplex/longest_repeat/stems/loops are recomputed here for every
    record from the actual sequence string, per the spec's step 3. Its
    `kd_nM` field IS reused as-is (spec marks it "already parsed").
  - AptaNexus's `Aptamer sequence` field is messy: ~4.6% of records carry
    oligo-synthesis annotations baked into the sequence text (5'/3' prime
    markers, fluorophore/quencher tags like FAM/BIOTIN/CY5, spacer codes
    like C3/C6/(CH2)6, LNA "+base" and phosphorothioate "*" linkage
    markers, per-base 2'-modification parentheticals). These are cleaned
    on a best-effort basis; sequences that still fail ATGCU-only
    validation afterward are dropped (counted and reported).
  - AptaNexus's `Affinity` field is free text (single value, multi-method
    comma-separated lists, ± error margins, ranges, or unparseable prose).
    Parsed per the spec's rules; falls back to `Kd = 10^(-pKd) * 1e9` when
    Affinity doesn't parse but pKd is present.
"""

import argparse
import json
import math
import re
from collections import Counter
from pathlib import Path

try:
    import RNA
    HAVE_VIENNARNA = True
except ImportError:
    HAVE_VIENNARNA = False


# ---------------------------------------------------------------------------
# Sequence cleaning
# ---------------------------------------------------------------------------

CYRILLIC_TO_LATIN = str.maketrans({
    "А": "A", "В": "B", "Е": "E", "К": "K", "М": "M", "Н": "H",
    "О": "O", "Р": "P", "С": "C", "Т": "T", "У": "Y", "Х": "X",
})

QUOTE_CLASS = r"['′‘’´]"
PRIME_MARKERS = re.compile(rf"[53]{QUOTE_CLASS}|{QUOTE_CLASS}[53]")
PAREN_ANNOTATION = re.compile(r"\([^)]*\)")

# Known oligo-synthesis modification tags (fluorophores, quenchers, linkers,
# spacers, end-modifiers). Every alternative here contains at least one
# character that is never valid in a pure ATGCU sequence, so stripping them
# can never remove content from an already-clean sequence.
MOD_TAGS = re.compile(
    r"BIOTIN|BIO|DITHIOL|THIOL|\bHS\b|"
    r"FAM|HEX|TET|ROX|JOE|CY5\.5|CY5|CY3|CY7|TAMRA|TEX615|ALEXA\d*|"
    r"DY\d{3}|IRDYE\d*|"
    r"BHQ\d|DABCYL|IABKFQ|IBFQ|"
    r"NH2|AMINO|PHOS(?:PHATE)?|"
    r"SPACER|SPC|PEG|TEG|HEG|MB|"
    r"\(CH2\)\d+|"
    r"C3|C6|C7|C9|C12|C15|C18",
    re.I,
)

# Whitespace/punctuation noise, plus LNA "+base" and phosphorothioate
# "*" linkage markers — stripping the marker character alone (not the base
# that follows) preserves the real nucleotide at that position.
STRIP_CHARS = re.compile(r"[\s\-\[\]:'\"´‘’′\+\*]")
EDGE_DIGIT = re.compile(r"^\d|\d$")
VALID_SEQ = re.compile(r"^[ATGCU]+$")


def clean_sequence(raw) -> str | None:
    if not isinstance(raw, str) or not raw.strip():
        return None
    s = raw.upper().translate(CYRILLIC_TO_LATIN)
    s = MOD_TAGS.sub("", s)
    s = PRIME_MARKERS.sub("", s)
    s = PAREN_ANNOTATION.sub("", s)
    s = STRIP_CHARS.sub("", s)
    while s and EDGE_DIGIT.search(s):
        s = s.strip("0123456789")
    if not s or not VALID_SEQ.match(s):
        return None
    return s


def detect_aptamer_type(sequence: str, fallback_label: str = "") -> str:
    has_u = "U" in sequence
    has_t = "T" in sequence
    if has_u and not has_t:
        return "RNA"
    if has_t and not has_u:
        return "DNA"
    return "RNA" if "RNA" in str(fallback_label).upper() else "DNA"


def clean_target_name(raw) -> str:
    if not isinstance(raw, str):
        return "Unknown"
    s = re.sub(r"\s+", " ", raw).strip().strip(".")
    return s if s else "Unknown"


# AptaDB's target_name is a raw UniProt entry code for ~47% of records
# (e.g. "THRB_HUMAN", "A6MI22_9HIV1") rather than a readable name — but a
# readable one is sitting right there in protein_name (e.g. "Prothrombin
# (EC 3.4.21.5) (Coagulation factor II) [Cleaved into: ...]"). Swap in the
# primary name (text before the first parenthetical/bracket) so the UI
# doesn't show cryptic codes and Target Lookup can group by a real name.
UNIPROT_CODE_PATTERN = re.compile(r"^[A-Z0-9]+_[A-Z0-9]+$")


def resolve_aptadb_target_name(target_name: str, protein_name) -> str:
    cleaned = clean_target_name(target_name)
    if UNIPROT_CODE_PATTERN.match(cleaned) and isinstance(protein_name, str) and protein_name.strip():
        primary = re.split(r"\s*[\(\[]", protein_name.strip(), maxsplit=1)[0].strip()
        if primary:
            return primary
    return cleaned


# Canonicalization for well-known targets that appear under many spelling
# variants across (and within) both sources — matters for Target Lookup,
# where fragmenting a famous target into near-duplicate entries would be
# misleading. Uses \b word boundaries so e.g. "Prothrombin" (a genuinely
# different protein) is deliberately left distinct from "Thrombin".
TARGET_ALIASES = [
    (re.compile(r"\bthrombin\b", re.I), "Thrombin"),
    (re.compile(r"vascular endothelial growth factor|\bvegf\b", re.I), "VEGF"),
    (re.compile(r"platelet-derived growth factor.*-?aa\b|pdgf-aa", re.I), "PDGF-AA"),
    (re.compile(r"platelet-derived growth factor.*-?bb\b|pdgf-bb", re.I), "PDGF-BB"),
    (re.compile(r"platelet-derived growth factor.*-?ab\b|pdgf-ab", re.I), "PDGF-AB"),
    (re.compile(r"reverse transcriptase.*hiv-1|hiv-1.*reverse transcriptase|hiv-1 rt", re.I), "HIV-1 Reverse Transcriptase"),
    (re.compile(r"rev protein of hiv-1|hiv-1 rev\b", re.I), "HIV-1 Rev"),
    (re.compile(r"taq\s*(dna)?\s*pol(ymerase)?|thermus aquaticus dna polymerase", re.I), "Taq DNA Polymerase"),
    (re.compile(r"immunoglobulin e|\bige\b", re.I), "IgE"),
    (re.compile(r"adenosine triphosphate|\batp\b", re.I), "ATP"),
    (re.compile(r"basic fibroblast growth factor|\bbfgf\b", re.I), "bFGF"),
]


def canonicalize_target(name: str) -> str:
    for pattern, canonical in TARGET_ALIASES:
        if pattern.search(name):
            return canonical
    return name


# ---------------------------------------------------------------------------
# Target type classification
# ---------------------------------------------------------------------------

def classify_aptanexus_target_type(raw: str) -> str:
    lower = (raw or "").lower()
    if "protein" in lower or "peptide" in lower:
        return "Protein"
    if "small molecule" in lower or "ion" in lower:
        return "Small Molecule"
    if "cell" in lower or "tissue" in lower:
        return "Cell"
    if "nucleic acid" in lower:
        return "Nucleic Acid"
    if "microorganism" in lower or "virus" in lower:
        return "Microorganism"
    return "Other"


APTADB_TARGET_TYPE_MAP = {
    "Protein": "Protein",
    "Molecule": "Small Molecule",
    "Cell": "Cell",
    "Other": "Other",
}


# ---------------------------------------------------------------------------
# Affinity / Kd parsing (AptaNexus's free-text "Affinity" field)
# ---------------------------------------------------------------------------

UNIT = r"(pM|nM|mM|u[Mm]|[μµ]M)"
KD_ERROR_MARGIN_PATTERN = re.compile(
    rf"([\d.]+)\s*(?:±|\+/-?|\+-)\s*[\d.]+\s*{UNIT}", re.I
)
KD_RANGE_PATTERN = re.compile(rf"([\d.]+)\s*[-–]\s*([\d.]+)\s*{UNIT}", re.I)
KD_PLAIN_PATTERN = re.compile(rf"([\d.]+)\s*{UNIT}", re.I)
KD_UNIT_TO_NM = {"pm": 1e-3, "nm": 1.0, "um": 1e3, "μm": 1e3, "µm": 1e3, "mm": 1e6}


# Real reported Kd values realistically span ~1 fM to ~100 mM. Values outside
# this band are near-certainly upstream unit/notation errors (e.g. AptaNexus
# rows whose Affinity text mixes scientific notation with a literal unit
# suffix, like "1.56 ×10⁻⁹ to 8.84 ×10⁻⁹ nM" — self-contradictory, and the
# matching pKd for the same row is corrupted too) rather than real
# measurements — reject rather than silently keep a number that misleads
# every downstream percentile/comparison.
KD_PLAUSIBLE_MIN_NM = 1e-3
KD_PLAUSIBLE_MAX_NM = 1e8


def _unit_to_nm(value: float, unit: str) -> float:
    return value * KD_UNIT_TO_NM[unit.lower()]


def _plausible_kd(value: float | None) -> float | None:
    if value is None:
        return None
    return value if KD_PLAUSIBLE_MIN_NM <= value <= KD_PLAUSIBLE_MAX_NM else None


def parse_kd_fragment(fragment: str) -> float | None:
    m = KD_ERROR_MARGIN_PATTERN.search(fragment)
    if m:
        return _unit_to_nm(float(m.group(1)), m.group(2))
    m = KD_RANGE_PATTERN.search(fragment)
    if m:
        lo, hi = float(m.group(1)), float(m.group(2))
        return _unit_to_nm((lo + hi) / 2, m.group(3))
    m = KD_PLAIN_PATTERN.search(fragment)
    if m:
        return _unit_to_nm(float(m.group(1)), m.group(2))
    return None


def parse_affinity_text(raw) -> float | None:
    if not isinstance(raw, str):
        return None
    s = raw.strip()
    if not s or s.upper() in ("N/A", "NA", "NOT PROVIDED", ""):
        return None
    values = [parse_kd_fragment(frag) for frag in s.split(",")]
    values = [_plausible_kd(v) for v in values if v is not None and v > 0]
    values = [v for v in values if v is not None]
    return min(values) if values else None


def kd_from_pkd(pkd) -> float | None:
    try:
        pkd = float(pkd)
    except (TypeError, ValueError):
        return None
    return _plausible_kd((10 ** (-pkd)) * 1e9)


def parse_year(raw) -> int | None:
    try:
        return int(float(raw))
    except (TypeError, ValueError):
        return None


def parse_doi(raw) -> str | None:
    if not isinstance(raw, str) or not raw.strip():
        return None
    return re.sub(r"^https?://(dx\.)?doi\.org/", "", raw.strip())


# ---------------------------------------------------------------------------
# Biophysical feature computation (applied uniformly to every record)
# ---------------------------------------------------------------------------

def dinucleotide_entropy(seq: str) -> float:
    if len(seq) < 2:
        return 0.0
    counts = Counter(seq[i:i + 2] for i in range(len(seq) - 1))
    total = sum(counts.values())
    probs = [c / total for c in counts.values()]
    return -sum(p * math.log2(p) for p in probs if p > 0)


G_QUAD_PATTERN = re.compile(r"G{3,}.{1,7}G{3,}.{1,7}G{3,}.{1,7}G{3,}")


def longest_homopolymer_run(seq: str) -> int:
    if not seq:
        return 0
    longest = current = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            current += 1
            longest = max(longest, current)
        else:
            current = 1
    return longest


def count_stems_loops(dot_bracket) -> tuple:
    if not isinstance(dot_bracket, str) or not dot_bracket:
        return None, None
    stems = dot_bracket.count("(")
    loops = len(re.findall(r"\(\.+\)", dot_bracket))
    return stems, loops


def predict_structure_if_missing(sequence: str, existing_mfe, existing_structure):
    """AptaDB records have no precomputed structure; fold them ourselves
    (same approach as v1) so the whole dataset has structure data where
    ViennaRNA is available. AptaNexus records already carry MFE/dot-bracket
    from the source and are left as-is."""
    if existing_mfe is not None or existing_structure:
        return existing_mfe, existing_structure
    if not HAVE_VIENNARNA or len(sequence) < 3:
        return None, None
    structure, mfe = RNA.fold(sequence.replace("T", "U"))
    return round(float(mfe), 2), structure


def compute_features(sequence: str, aptamer_type: str, mfe_structure) -> dict:
    length = len(sequence)
    counts = Counter(sequence)
    a, g, c = counts.get("A", 0), counts.get("G", 0), counts.get("C", 0)
    t_or_u = counts.get("U" if aptamer_type == "RNA" else "T", 0)
    stems, loops = count_stems_loops(mfe_structure)

    return {
        "gc_content": round((g + c) / length, 4),
        "a_freq": round(a / length, 4),
        "t_freq": round(t_or_u / length, 4),
        "g_freq": round(g / length, 4),
        "c_freq": round(c / length, 4),
        "purine_ratio": round((a + g) / length, 4),
        "complexity": round(dinucleotide_entropy(sequence), 4),
        "has_g_quadruplex": bool(G_QUAD_PATTERN.search(sequence)),
        "longest_repeat": longest_homopolymer_run(sequence),
        "num_stems": stems,
        "num_loops": loops,
    }


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

def load_aptanexus(path: Path) -> list[dict]:
    records = []
    dropped_no_seq = 0
    dropped_invalid = 0

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            row = json.loads(line)

            raw_seq = row.get("Aptamer sequence")
            if not isinstance(raw_seq, str) or not raw_seq.strip():
                dropped_no_seq += 1
                continue

            sequence = clean_sequence(raw_seq)
            if sequence is None:
                dropped_invalid += 1
                continue

            aptamer_type = detect_aptamer_type(sequence)
            target_name = canonicalize_target(clean_target_name(row.get("Target name")))
            target_type = classify_aptanexus_target_type(row.get("Target type"))

            kd_nM = parse_affinity_text(row.get("Affinity"))
            pkd = row.get("pKd")
            pkd = float(pkd) if isinstance(pkd, (int, float)) else None
            if kd_nM is None and pkd is not None:
                kd_nM = kd_from_pkd(pkd)

            mfe_structure = row.get("SecStr_dotbracket") or None
            predicted_mfe = row.get("MFE")
            predicted_mfe = float(predicted_mfe) if isinstance(predicted_mfe, (int, float)) else None

            record = {
                "sequence": sequence,
                "aptamer_type": aptamer_type,
                "length": len(sequence),
                "target_name": target_name,
                "target_type": target_type,
                "kd_nM": kd_nM,
                "pkd": pkd,
                "predicted_mfe": predicted_mfe,
                "mfe_structure": mfe_structure,
                "quality_tier": row.get("Level") or None,
                "year": parse_year(row.get("Year")),
                "doi": parse_doi(row.get("Doi")),
                "journal": row.get("Journal") or None,
                "article_title": row.get("Article title") or None,
                "external_id": row.get("External_ID") or None,
                "external_name": row.get("External_Name") or None,
                "gene_symbol": row.get("Gene_Symbol") or None,
                "id_type": row.get("ID_Type") or None,
                "buffer": row.get("Buffer condition") or None,
                "is_best": row.get("Best") if isinstance(row.get("Best"), bool) else None,
                "sequence_id": row.get("Sequence ID") or None,
                "uniprot_id": None,
                "protein_name": None,
                "pubmed_id": None,
                "selex_method": None,
                "source": "aptanexus",
            }
            records.append(record)

    print(f"[aptanexus] Loaded rows: {dropped_no_seq + dropped_invalid + len(records)}")
    print(f"[aptanexus] Dropped (no sequence): {dropped_no_seq}")
    print(f"[aptanexus] Dropped (invalid chars after cleaning): {dropped_invalid}")
    print(f"[aptanexus] Remaining: {len(records)}")
    return records


def load_aptadb(path: Path) -> list[dict]:
    with open(path) as f:
        rows = json.load(f)

    records = []
    dropped_no_seq = 0
    dropped_invalid = 0

    for row in rows:
        raw_seq = row.get("sequence")
        if not isinstance(raw_seq, str) or not raw_seq.strip():
            dropped_no_seq += 1
            continue

        sequence = clean_sequence(raw_seq)
        if sequence is None:
            dropped_invalid += 1
            continue

        aptamer_type = detect_aptamer_type(sequence, row.get("aptamer_type", ""))
        target_name = canonicalize_target(
            resolve_aptadb_target_name(row.get("target_name"), row.get("protein_name"))
        )
        target_type = APTADB_TARGET_TYPE_MAP.get(row.get("target_chemistry"), "Other")

        kd_nM = row.get("kd_nM")
        kd_nM = _plausible_kd(float(kd_nM)) if isinstance(kd_nM, (int, float)) and kd_nM > 0 else None

        record = {
            "sequence": sequence,
            "aptamer_type": aptamer_type,
            "length": len(sequence),
            "target_name": target_name,
            "target_type": target_type,
            "kd_nM": kd_nM,
            "pkd": None,
            "predicted_mfe": None,  # filled in later if ViennaRNA available
            "mfe_structure": row.get("mfe_structure") or None,
            "quality_tier": None,
            "year": None,
            "doi": None,
            "journal": None,
            "article_title": None,
            "external_id": None,
            "external_name": None,
            "gene_symbol": row.get("gene_names") or None,
            "id_type": "UniProt ID" if row.get("uniprot_id") else None,
            "buffer": row.get("buffer_conditions") or None,
            "is_best": None,
            "sequence_id": row.get("aptamer_id") or None,
            "uniprot_id": row.get("uniprot_id") or None,
            "protein_name": row.get("protein_name") or None,
            "pubmed_id": row.get("pubmed_id") or None,
            "selex_method": None,
            "source": "aptadb",
        }
        records.append(record)

    print(f"[aptadb] Loaded rows: {dropped_no_seq + dropped_invalid + len(records)}")
    print(f"[aptadb] Dropped (no sequence): {dropped_no_seq}")
    print(f"[aptadb] Dropped (invalid chars after cleaning): {dropped_invalid}")
    print(f"[aptadb] Remaining: {len(records)}")
    return records


# ---------------------------------------------------------------------------
# Dedup + enrich + output
# ---------------------------------------------------------------------------

def completeness_score(r: dict) -> int:
    return sum(1 for v in r.values() if v is not None and v != "")


def deduplicate(records: list[dict]) -> list[dict]:
    best_by_key = {}
    dup_keys = set()
    for r in records:
        key = (r["sequence"], r["target_name"])
        existing = best_by_key.get(key)
        if existing is None:
            best_by_key[key] = r
        else:
            dup_keys.add(key)
            winner = r if completeness_score(r) > completeness_score(existing) else existing
            best_by_key[key] = winner

    deduped = list(best_by_key.values())
    for r in deduped:
        key = (r["sequence"], r["target_name"])
        if key in dup_keys:
            r["source"] = "both"

    print(f"Deduplicated {len(records)} -> {len(deduped)} records ({len(dup_keys)} overlapping keys)")
    return deduped


def enrich(records: list[dict]) -> list[dict]:
    enriched = []
    for i, r in enumerate(records):
        predicted_mfe, mfe_structure = predict_structure_if_missing(
            r["sequence"], r["predicted_mfe"], r["mfe_structure"]
        )
        r["predicted_mfe"] = predicted_mfe
        r["mfe_structure"] = mfe_structure
        features = compute_features(r["sequence"], r["aptamer_type"], mfe_structure)
        full_record = {"id": f"{r['source']}_{i + 1:05d}" if r["source"] != "both" else f"apta_{i + 1:05d}",
                        **r, **features}
        enriched.append(full_record)
    return enriched


def compute_stats(records: list[dict]) -> dict:
    kds = [r["kd_nM"] for r in records if r["kd_nM"] is not None]
    lengths = sorted(r["length"] for r in records)
    gcs = sorted(r["gc_content"] for r in records)

    def median(vals):
        n = len(vals)
        if n == 0:
            return None
        mid = n // 2
        return vals[mid] if n % 2 else (vals[mid - 1] + vals[mid]) / 2

    return {
        "total_records": len(records),
        "records_with_kd": len(kds),
        "records_with_structure": sum(1 for r in records if r["mfe_structure"]),
        "unique_targets": len({r["target_name"] for r in records}),
        "dna_count": sum(1 for r in records if r["aptamer_type"] == "DNA"),
        "rna_count": sum(1 for r in records if r["aptamer_type"] == "RNA"),
        "kd_range": [min(kds), max(kds)] if kds else [None, None],
        "median_length": median(lengths),
        "median_gc": median(gcs),
        "source_counts": dict(Counter(r["source"] for r in records)),
        "quality_tier_counts": dict(Counter(r["quality_tier"] for r in records if r["quality_tier"])),
        "vienna_rna_available": HAVE_VIENNARNA,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--aptanexus", required=True)
    parser.add_argument("--aptadb", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--stats-output", default=None)
    args = parser.parse_args()

    records = load_aptanexus(Path(args.aptanexus)) + load_aptadb(Path(args.aptadb))
    records = deduplicate(records)
    records = enrich(records)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(records, f, separators=(',', ':'))
    print(f"Wrote {len(records)} records to {out_path}")

    stats = compute_stats(records)
    stats_path = Path(args.stats_output) if args.stats_output else out_path.parent / "dataset_stats.json"
    with open(stats_path, "w") as f:
        json.dump(stats, f, indent=2)
    print(f"Wrote stats to {stats_path}")
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
