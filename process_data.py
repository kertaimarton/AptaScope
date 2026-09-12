"""
AptaScope — Phase 1 data processing.

Reads two real aptamer datasets:
  - UTexas Aptamer Database (Zenodo record 8387047, "UTexas Aptamer Database
    dataset_Sept2023.xlsx") — 1,495 rows with sequence, target, Kd, DOI, year.
  - Ribocentre Aptamer (aptamer.ribocentre.org) — 509 rows, fetched from its
    underlying data endpoint (apidata/sequences_cleaned.json), with sequence,
    target ("Ligand"), a direct target-type label, free-text affinity, year.

and produces:
  - src/data/aptamer_data.json   (array of enriched aptamer records)
  - src/data/dataset_stats.json  (summary stats)

Notes on sourcing:
AptaDB (lmmd.ecust.edu.cn/aptadb) renders its search results client-side via
JavaScript with no discoverable static/JSON endpoint reachable by a plain
HTTP client, so it could not be scraped without a headless browser. Per the
spec's fallback instructions, this script substitutes Ribocentre Aptamer
instead (also a curated literature database, directly downloadable), on top
of the UTexas priority source. Both are real, richly annotated experimental
data — not mock data. The two sources overlap substantially (much of the
same SELEX literature is catalogued in both); overlapping records are
resolved by the dedup step in favor of whichever has the more complete
annotation (UTexas rows carry a real DOI, which usually wins the tie-break).
"""

import json
import math
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

try:
    import RNA
    HAVE_VIENNARNA = True
except ImportError:
    HAVE_VIENNARNA = False

RAW_XLSX = Path(__file__).parent / "raw_data" / "UTexas_Aptamer_Database_Sept2023.xlsx"
RAW_RIBOCENTRE = Path(__file__).parent / "raw_data" / "ribocentre_sequences.json"
OUT_DIR = Path(__file__).parent / "src" / "data"
OUT_DATA = OUT_DIR / "aptamer_data.json"
OUT_STATS = OUT_DIR / "dataset_stats.json"

# ---------------------------------------------------------------------------
# Step 1.2 — Clean and unify
# ---------------------------------------------------------------------------

# Homoglyph characters that show up from copy-paste errors in the source
# spreadsheet (Cyrillic look-alikes mixed into otherwise-Latin sequences).
CYRILLIC_TO_LATIN = str.maketrans({
    "А": "A", "В": "B", "Е": "E", "К": "K", "М": "M", "Н": "H",
    "О": "O", "Р": "P", "С": "C", "Т": "T", "У": "Y", "Х": "X",
})

QUOTE_CLASS = r"['′‘’´]"
# Matches 5'/3' markers in either digit-then-quote or quote-then-digit order
# (a handful of source rows write the terminus marker reversed, e.g. "...ACGT'3"),
# and accepts the acute-accent character some rows use in place of a prime.
PRIME_MARKERS = re.compile(rf"[53]{QUOTE_CLASS}|{QUOTE_CLASS}[53]")
PAREN_ANNOTATION = re.compile(r"\([^)]*\)")
STRIP_CHARS = re.compile(r"[\s\-\[\]:'\"´‘’′]")
VALID_SEQ = re.compile(r"^[ATGCU]+$")
EDGE_DIGIT = re.compile(r"^\d|\d$")


def clean_sequence(raw: str) -> str | None:
    if not isinstance(raw, str):
        return None
    s = raw.upper().translate(CYRILLIC_TO_LATIN)
    s = PRIME_MARKERS.sub("", s)
    s = PAREN_ANNOTATION.sub("", s)  # drop modification annotations e.g. (2'FC)
    s = STRIP_CHARS.sub("", s)
    # Safety net: a stray leading/trailing digit left over from a terminus
    # marker whose quote character got stripped without being recognized
    # above (e.g. a bare "...SEQ3" with no accompanying quote at all).
    while s and EDGE_DIGIT.search(s):
        s = s.strip("0123456789")
    if not s or not VALID_SEQ.match(s):
        return None
    return s


def detect_type(sequence: str, fallback_label: str) -> str:
    has_u = "U" in sequence
    has_t = "T" in sequence
    if has_u and not has_t:
        return "RNA"
    if has_t and not has_u:
        return "DNA"
    return "RNA" if "RNA" in str(fallback_label).upper() else "DNA"


def clean_target_name(raw: str) -> str:
    if not isinstance(raw, str):
        return "Unknown"
    s = re.sub(r"\s+", " ", raw).strip().strip(".")
    return s if s else "Unknown"


# Canonicalization for well-known targets that appear under many spelling
# variants in the source data (matters for the Target Lookup view, where
# fragmenting "Thrombin, Human" / "alpha-Thrombin, Human" / etc. into
# separate targets would be misleading).
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


SMALL_MOLECULE_KEYWORDS = [
    "atp", "adenosine", "cocaine", "theophylline", "dopamine", "serotonin",
    "chloramphenicol", "kanamycin", "streptomycin", "tobramycin", "neomycin",
    "lividomycin", "ochratoxin", "aflatoxin", "melamine", "sulfadimethoxine",
    "bisphenol", "estradiol", "testosterone", "cortisol", "vitamin", "biotin",
    "folic acid", "glucose", "sialyllactose", "cholesterol", "cyanocobalamin",
    "cobinamide", "citrulline", "arginine", "sphingosylphosphorylcholine",
    "cibacron blue", "reactive blue", "hematoporphyrin", "nmm", "porphyrin",
    "mesoporphyrin", "malachite green", "rhodamine", "fluorescein", "caffeine",
    "histamine", "epinephrine", "norepinephrine", "cortisone", "progesterone",
    "digoxin", "warfarin", "penicillin", "ampicillin", "tetracycline",
    "sulfonamide", "phthalate", "toxin b", "microcystin", "atrazine",
]

CELL_KEYWORDS = [
    "cell line", "cells,", " cells", "hepg2", "ccrf-cem", "ramos", "a549",
    "hela", "jurkat", "mcf-7", "hpac", "pl45", "bnl 1me", "lymphocyte",
    "carcinoma cells", "cancer cell", "tumor cell",
]

ORGANISM_KEYWORDS = [
    "virus", "bacteri", "coli", "salmonella", "streptococcus", "staphylococcus",
    "parasite", "oocyst", "spore", "fungus", "cryptosporidium", "mycobacterium",
    "listeria", "iridovirus", "prion",
]

PROTEIN_KEYWORDS = [
    "protein", "enzyme", "kinase", "polymerase", "receptor", "antigen",
    "antibody", "factor", "immunoglobulin", "albumin", "lysozyme",
    "interferon", "interleukin", "synthetase", "integrase", "transcriptase",
    "peptide", "hormone", "growth factor", "cytokine", "nuclease", "protease",
    "phosphatase", "reductase", "channel", "transporter", "collagen",
]


def parse_molecular_weight(raw) -> float | None:
    if not isinstance(raw, str):
        return None
    m = re.search(r"([\d.]+)\s*Da", raw, re.I)
    return float(m.group(1)) if m else None


def classify_target_type(name: str, mw_raw) -> str:
    lower = name.lower()
    mw = parse_molecular_weight(mw_raw)

    if any(k in lower for k in CELL_KEYWORDS):
        return "cell"
    if any(k in lower for k in ORGANISM_KEYWORDS):
        return "other"
    if any(k in lower for k in SMALL_MOLECULE_KEYWORDS):
        return "small_molecule"
    if mw is not None and mw < 900:
        return "small_molecule"
    if any(k in lower for k in PROTEIN_KEYWORDS):
        return "protein"
    if mw is not None and mw >= 5000:
        return "protein"
    return "protein"  # majority-class fallback: most SELEX targets are proteins


def parse_kd(raw) -> float | None:
    val = pd.to_numeric(raw, errors="coerce")
    if pd.isna(val) or val <= 0:
        return None
    return float(val)


def parse_year(raw) -> int | None:
    val = pd.to_numeric(raw, errors="coerce")
    return int(val) if not pd.isna(val) else None


def parse_doi(raw) -> str | None:
    if not isinstance(raw, str):
        return None
    s = raw.strip()
    if not s or s.lower() in ("not reported", "n/a", "none"):
        return None
    return re.sub(r"^https?://(dx\.)?doi\.org/", "", s)


# Ribocentre reports Kd as free text (e.g. "Kd: 23±3 nM", "Kd: <100µM",
# "Kd1: 0.19±0.02 nM\nKd2: 49±26 nM"). Extract the first clean number+unit
# pair, tolerating an inline "±error" between the value and its unit; report
# text with no parseable number/unit (rate constants in 1/s, percentages,
# scientific-notation Molar, "NA") is left as null rather than guessed at.
KD_TEXT_PATTERN = re.compile(r"([\d.]+)(?:\s*±\s*[\d.]+)?\s*(pM|nM|[μµ]M|mM)", re.I)
KD_UNIT_TO_NM = {"pm": 1e-3, "nm": 1.0, "μm": 1e3, "µm": 1e3, "mm": 1e6}


def parse_kd_text(raw) -> float | None:
    if not isinstance(raw, str):
        return None
    m = KD_TEXT_PATTERN.search(raw)
    if not m:
        return None
    value, unit = float(m.group(1)), m.group(2).lower()
    kd_nM = value * KD_UNIT_TO_NM[unit]
    return kd_nM if kd_nM > 0 else None


# ---------------------------------------------------------------------------
# Step 1.3 — Biophysical features
# ---------------------------------------------------------------------------

def dinucleotide_entropy(seq: str) -> float:
    if len(seq) < 2:
        return 0.0
    dinucs = [seq[i:i + 2] for i in range(len(seq) - 1)]
    counts = Counter(dinucs)
    total = sum(counts.values())
    probs = [c / total for c in counts.values()]
    return -sum(p * math.log2(p) for p in probs if p > 0)


G_QUAD_PATTERN = re.compile(r"G{3,}.{1,7}G{3,}.{1,7}G{3,}.{1,7}G{3,}")


def has_g_quadruplex(seq: str) -> bool:
    return bool(G_QUAD_PATTERN.search(seq))


def longest_homopolymer_run(seq: str) -> int:
    longest = 1
    current = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            current += 1
            longest = max(longest, current)
        else:
            current = 1
    return longest if seq else 0


def predict_structure(seq: str):
    if not HAVE_VIENNARNA or len(seq) < 3:
        return None, None, None
    fold_seq = seq.replace("T", "U")
    structure, mfe = RNA.fold(fold_seq)
    stems = structure.count("(")
    loops = len(re.findall(r"\(\.+\)", structure))
    return round(float(mfe), 2), stems, loops


def compute_features(sequence: str, seq_type: str) -> dict:
    length = len(sequence)
    counts = Counter(sequence)
    a = counts.get("A", 0)
    t_or_u = counts.get("U" if seq_type == "RNA" else "T", 0)
    g = counts.get("G", 0)
    c = counts.get("C", 0)

    mfe, stems, loops = predict_structure(sequence)

    return {
        "gc_content": round((g + c) / length, 4),
        "a_freq": round(a / length, 4),
        "t_freq": round(t_or_u / length, 4),
        "g_freq": round(g / length, 4),
        "c_freq": round(c / length, 4),
        "purine_ratio": round((a + g) / length, 4),
        "complexity": round(dinucleotide_entropy(sequence), 4),
        "has_g_quadruplex": has_g_quadruplex(sequence),
        "longest_repeat": longest_homopolymer_run(sequence),
        "predicted_mfe": mfe,
        "num_stems": stems,
        "num_loops": loops,
    }


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def load_utexas() -> list[dict]:
    df = pd.read_excel(RAW_XLSX, sheet_name="Dataset")

    records = []
    dropped_no_sequence = 0
    dropped_invalid_chars = 0

    for _, row in df.iterrows():
        raw_seq = row["Aptamer Sequence"]
        if not isinstance(raw_seq, str) or not raw_seq.strip():
            dropped_no_sequence += 1
            continue

        sequence = clean_sequence(raw_seq)
        if sequence is None:
            dropped_invalid_chars += 1
            continue

        seq_type = detect_type(sequence, row.get("Type of Nucleic Acid", ""))
        target_name = canonicalize_target(clean_target_name(row.get("Target ")))

        record = {
            "sequence": sequence,
            "type": seq_type,
            "length": len(sequence),
            "target_name": target_name,
            "target_type": classify_target_type(target_name, row.get("Molecular weight of target")),
            "kd_nM": parse_kd(row.get("Kd (nM)")),
            "selex_method": None,
            "year": parse_year(row.get("Year of Paper")),
            "doi": parse_doi(row.get("Journal DOI")),
            "source": "utexas",
        }
        records.append(record)

    print(f"[utexas] Loaded {len(df)} raw rows")
    print(f"[utexas] Dropped (no sequence): {dropped_no_sequence}")
    print(f"[utexas] Dropped (invalid characters after cleaning): {dropped_invalid_chars}")
    print(f"[utexas] Remaining: {len(records)}")

    return records


# Direct target-type label provided by Ribocentre; more reliable than the
# keyword heuristic (which exists because UTexas doesn't provide this field).
RIBOCENTRE_TYPE_MAP = {
    "Proteins": "protein",
    "Small molecules": "small_molecule",
    "Cells": "cell",
    "Nucleic acids": "other",
}


def load_ribocentre() -> list[dict]:
    if not RAW_RIBOCENTRE.exists():
        print("[ribocentre] data file not found, skipping")
        return []

    with open(RAW_RIBOCENTRE) as f:
        rows = json.load(f)["Sheet1"]

    records = []
    dropped_no_sequence = 0
    dropped_invalid_chars = 0

    for row in rows:
        raw_seq = row.get("Sequence")
        if not isinstance(raw_seq, str) or not raw_seq.strip():
            dropped_no_sequence += 1
            continue

        sequence = clean_sequence(raw_seq)
        if sequence is None:
            dropped_invalid_chars += 1
            continue

        seq_type = detect_type(sequence, "")
        target_name = canonicalize_target(clean_target_name(row.get("Ligand")))
        target_type = RIBOCENTRE_TYPE_MAP.get(row.get("Type"), "other")

        record = {
            "sequence": sequence,
            "type": seq_type,
            "length": len(sequence),
            "target_name": target_name,
            "target_type": target_type,
            "kd_nM": parse_kd_text(row.get("Affinity")),
            "selex_method": None,
            "year": parse_year(row.get("Year")),
            "doi": None,
            "source": "ribocentre",
        }
        records.append(record)

    print(f"[ribocentre] Loaded {len(rows)} raw rows")
    print(f"[ribocentre] Dropped (no sequence): {dropped_no_sequence}")
    print(f"[ribocentre] Dropped (invalid characters after cleaning): {dropped_invalid_chars}")
    print(f"[ribocentre] Remaining: {len(records)}")

    return records


def completeness_score(r: dict) -> tuple:
    return (
        r["kd_nM"] is not None,
        r["doi"] is not None,
        r["year"] is not None,
        len(r["target_name"]),
    )


def deduplicate(records: list[dict]) -> list[dict]:
    best_by_key = {}
    for r in records:
        key = (r["sequence"], r["target_name"])
        existing = best_by_key.get(key)
        if existing is None or completeness_score(r) > completeness_score(existing):
            best_by_key[key] = r
    deduped = list(best_by_key.values())
    print(f"Deduplicated {len(records)} -> {len(deduped)} records")
    return deduped


def enrich(records: list[dict]) -> list[dict]:
    enriched = []
    for i, r in enumerate(records):
        features = compute_features(r["sequence"], r["type"])
        full_record = {"id": f"APTA-{i + 1:05d}", **r, **features}
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
        "unique_targets": len({r["target_name"] for r in records}),
        "dna_count": sum(1 for r in records if r["type"] == "DNA"),
        "rna_count": sum(1 for r in records if r["type"] == "RNA"),
        "kd_range": [min(kds), max(kds)] if kds else [None, None],
        "median_length": median(lengths),
        "median_gc": median(gcs),
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "vienna_rna_available": HAVE_VIENNARNA,
    }


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    records = load_utexas() + load_ribocentre()
    records = deduplicate(records)
    records = enrich(records)

    with open(OUT_DATA, "w") as f:
        json.dump(records, f, indent=2)
    print(f"Wrote {len(records)} records to {OUT_DATA}")

    stats = compute_stats(records)
    with open(OUT_STATS, "w") as f:
        json.dump(stats, f, indent=2)
    print(f"Wrote stats to {OUT_STATS}")
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
