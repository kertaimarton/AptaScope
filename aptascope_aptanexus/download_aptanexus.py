"""
AptaNexus Downloader — Downloads the full AptaNexus aptamer database
https://www.aptanexus.com/

AptaNexus contains 12,500+ curated records, 4,815 affinity-validated
sequences, 1,900+ unique targets, mined from 23,000+ publications (2005-2025).

The database is available as a single JSONL file (one JSON object per line).

Outputs:
  - aptanexus_raw.jsonl   (raw download, preserved as-is)
  - aptanexus_raw.json    (parsed into a JSON array)
  - aptanexus_raw.csv     (flat CSV for quick inspection)

Usage:
  pip install requests --break-system-packages
  python download_aptanexus.py
"""

import requests
import json
import csv
import os
import sys
import re

# Primary: bulk JSONL download (seen on AptaNexus homepage)
JSONL_URL = "https://aptamer-database.vercel.app/APTAMERS.jsonl"

# Fallback: try the main site's API patterns
SEARCH_API_CANDIDATES = [
    "https://www.aptanexus.com/api/aptamers",
    "https://www.aptanexus.com/api/search",
    "https://www.aptanexus.com/api/data",
    "https://aptamer-database.vercel.app/api/aptamers",
    "https://aptamer-database.vercel.app/api/search",
]

OUTPUT_JSONL = "aptanexus_raw.jsonl"
OUTPUT_JSON = "aptanexus_raw.json"
OUTPUT_CSV = "aptanexus_raw.csv"
TIMEOUT = 60


def download_jsonl(url: str) -> list[dict]:
    """Download and parse a JSONL file (one JSON object per line)."""
    print(f"  Downloading from {url} ...")

    response = requests.get(url, timeout=TIMEOUT, headers={
        "User-Agent": "AptaScope-Research/1.0 (academic research tool)",
        "Accept": "application/jsonl, application/json, text/plain, */*",
    })

    if response.status_code != 200:
        print(f"  HTTP {response.status_code} — download failed")
        return []

    content = response.text.strip()

    if not content:
        print("  Empty response")
        return []

    records = []
    errors = 0

    for i, line in enumerate(content.split("\n")):
        line = line.strip()
        if not line:
            continue
        try:
            record = json.loads(line)
            records.append(record)
        except json.JSONDecodeError:
            errors += 1
            if errors <= 3:
                print(f"  Line {i+1}: failed to parse JSON — {line[:80]}...")

    print(f"  Parsed {len(records)} records ({errors} parse errors)")
    return records


def try_json_api(url: str) -> list[dict]:
    """Try fetching from a JSON API endpoint."""
    print(f"  Trying {url} ...")
    try:
        response = requests.get(url, timeout=TIMEOUT, headers={
            "User-Agent": "AptaScope-Research/1.0",
            "Accept": "application/json",
        })

        if response.status_code != 200:
            print(f"  HTTP {response.status_code}")
            return []

        data = response.json()

        # Handle different response shapes
        if isinstance(data, list):
            print(f"  Got array of {len(data)} records")
            return data
        elif isinstance(data, dict):
            # Look for common data keys
            for key in ["data", "results", "records", "aptamers", "items"]:
                if key in data and isinstance(data[key], list):
                    print(f"  Got {len(data[key])} records under '{key}'")
                    return data[key]
            # Maybe paginated — check for total/count
            if "total" in data or "count" in data:
                total = data.get("total", data.get("count", "?"))
                print(f"  Paginated response — total: {total}")
                return fetch_paginated(url, data)

        print(f"  Unexpected response shape: {type(data)}")
        return []

    except requests.exceptions.ConnectionError:
        print(f"  Connection error")
        return []
    except requests.exceptions.Timeout:
        print(f"  Timeout")
        return []
    except json.JSONDecodeError:
        print(f"  Response is not JSON")
        return []


def fetch_paginated(base_url: str, first_response: dict) -> list[dict]:
    """Handle paginated API responses."""
    all_records = []

    # Extract records from first page
    for key in ["data", "results", "records", "aptamers", "items"]:
        if key in first_response and isinstance(first_response[key], list):
            all_records.extend(first_response[key])
            break

    total = first_response.get("total", first_response.get("count", 0))
    if not total or total <= len(all_records):
        return all_records

    # Determine pagination params
    page_size = len(all_records) if all_records else 100
    pages_needed = (total // page_size) + 1

    print(f"  Fetching {pages_needed} pages (total: {total}, page_size: {page_size})")

    for page in range(1, min(pages_needed + 1, 200)):  # cap at 200 pages
        # Try common pagination patterns
        for param_style in [
            f"?page={page+1}&limit={page_size}",
            f"?offset={page * page_size}&limit={page_size}",
            f"?page={page+1}&size={page_size}",
            f"?skip={page * page_size}&take={page_size}",
        ]:
            try:
                url = base_url + param_style
                resp = requests.get(url, timeout=TIMEOUT, headers={
                    "User-Agent": "AptaScope-Research/1.0",
                    "Accept": "application/json",
                })
                if resp.status_code == 200:
                    data = resp.json()
                    page_records = []
                    if isinstance(data, list):
                        page_records = data
                    elif isinstance(data, dict):
                        for key in ["data", "results", "records", "aptamers", "items"]:
                            if key in data and isinstance(data[key], list):
                                page_records = data[key]
                                break

                    if page_records:
                        all_records.extend(page_records)
                        print(f"  Page {page+1}: +{len(page_records)} records (total: {len(all_records)})")
                        break
                    elif len(page_records) == 0:
                        print(f"  Page {page+1}: empty — likely done")
                        return all_records
            except Exception:
                continue

        # Rate limit
        import time
        time.sleep(0.5)

    return all_records


def normalize_record(record: dict) -> dict:
    """
    Normalize an AptaNexus record into the common AptaScope schema.
    Since we don't know the exact field names yet, this handles
    multiple possible naming conventions.
    """
    normalized = {}

    # Map common field names to our schema
    field_maps = {
        "sequence": ["sequence", "aptamer_sequence", "seq", "Sequence", "aptamerSequence"],
        "target_name": ["target", "target_name", "targetName", "Target", "target_molecule"],
        "target_type": ["target_type", "targetType", "target_category", "category", "targetCategory"],
        "kd_raw": ["kd", "Kd", "KD", "affinity", "binding_affinity", "bindingAffinity", "dissociation_constant"],
        "aptamer_type": ["type", "aptamer_type", "nucleic_acid_type", "chemistry", "aptamerType", "nucleicAcidType"],
        "length": ["length", "seq_length", "sequence_length", "Length"],
        "gc_content": ["gc_content", "gc", "GC", "gcContent", "gc_percentage"],
        "pubmed_id": ["pubmed_id", "pmid", "PubMed", "pubmedId", "PMID", "pubmed"],
        "doi": ["doi", "DOI", "paper_doi"],
        "year": ["year", "Year", "publication_year", "pub_year"],
        "quality_tier": ["tier", "quality_tier", "quality", "qualityTier", "data_quality"],
        "selex_method": ["selex", "selex_method", "selexMethod", "selection_method"],
        "buffer": ["buffer", "binding_buffer", "bindingBuffer", "conditions"],
        "modification": ["modification", "modifications", "chemical_modification"],
        "application": ["application", "applications", "use", "function"],
    }

    for target_key, source_keys in field_maps.items():
        for sk in source_keys:
            if sk in record and record[sk] is not None and str(record[sk]).strip():
                normalized[target_key] = record[sk]
                break

    # Parse affinity to nM if present
    kd_raw = normalized.get("kd_raw")
    if kd_raw is not None:
        normalized["kd_nM"] = parse_affinity(kd_raw)

    # Compute length if sequence present but length missing
    if "sequence" in normalized and "length" not in normalized:
        clean_seq = re.sub(r"[^ATGCUatgcu]", "", str(normalized["sequence"]))
        normalized["length"] = len(clean_seq)

    # Detect DNA vs RNA if not specified
    if "sequence" in normalized and "aptamer_type" not in normalized:
        seq = str(normalized["sequence"]).upper()
        if "U" in seq and "T" not in seq:
            normalized["aptamer_type"] = "RNA"
        elif "T" in seq:
            normalized["aptamer_type"] = "DNA"

    # Clean sequence
    if "sequence" in normalized:
        normalized["sequence"] = re.sub(r"[^ATGCUatgcu]", "", str(normalized["sequence"]).upper())

    # Keep all original fields too (prefixed) for inspection
    normalized["_raw"] = record

    normalized["source"] = "aptanexus"

    return normalized


def parse_affinity(value) -> float | None:
    """Convert various affinity representations to nanomolar."""
    if isinstance(value, (int, float)):
        # If it's a bare number, assume nM (most common in databases)
        return float(value)

    if not isinstance(value, str):
        return None

    value = value.strip()

    match = re.search(r"([\d.]+)\s*(nM|μM|uM|mM|pM|fM|M|nm|um|mm|pm|fm)", value, re.IGNORECASE)
    if not match:
        # Try bare number
        try:
            return float(value)
        except ValueError:
            return None

    num = float(match.group(1))
    unit = match.group(2).lower()

    multipliers = {
        "fm": 1e-6, "pm": 1e-3, "nm": 1.0, "μm": 1e3,
        "um": 1e3, "mm": 1e6, "m": 1e9,
    }

    return round(num * multipliers.get(unit, 1.0), 4)


def save_outputs(records: list[dict], raw_records: list[dict]):
    """Save processed and raw outputs."""

    # Save raw JSONL (preserve original)
    with open(OUTPUT_JSONL, "w", encoding="utf-8") as f:
        for r in raw_records:
            f.write(json.dumps(r, ensure_ascii=False) + "\n")
    print(f"  Saved raw JSONL: {OUTPUT_JSONL}")

    # Save normalized JSON
    with open(OUTPUT_JSON, "w", encoding="utf-8") as f:
        # Strip the _raw field for the clean version
        clean = []
        for r in records:
            c = {k: v for k, v in r.items() if k != "_raw"}
            clean.append(c)
        json.dump(clean, f, indent=2, ensure_ascii=False)
    print(f"  Saved normalized JSON: {OUTPUT_JSON}")

    # Save CSV (flat fields only)
    if not records:
        return

    flat_keys = set()
    for r in records:
        for k, v in r.items():
            if k != "_raw" and not isinstance(v, (list, dict)):
                flat_keys.add(k)
    flat_keys = sorted(flat_keys)

    with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=flat_keys, extrasaction="ignore")
        writer.writeheader()
        for r in records:
            writer.writerow({k: r.get(k, "") for k in flat_keys})
    print(f"  Saved CSV: {OUTPUT_CSV}")


def main():
    print("=" * 60)
    print("AptaNexus Downloader for AptaScope")
    print("Source: https://www.aptanexus.com/")
    print("Expected: ~12,500 aptamer records")
    print("=" * 60)

    raw_records = []

    # Strategy 1: Direct JSONL download (best case)
    print(f"\n[1/3] Trying bulk JSONL download...")
    raw_records = download_jsonl(JSONL_URL)

    # Strategy 2: Try API endpoints
    if not raw_records:
        print(f"\n[2/3] JSONL failed — trying API endpoints...")
        for url in SEARCH_API_CANDIDATES:
            raw_records = try_json_api(url)
            if raw_records:
                break

    # Strategy 3: Try search with wildcard/empty query
    if not raw_records:
        print(f"\n[3/3] APIs failed — trying search with broad queries...")
        search_urls = [
            "https://www.aptanexus.com/api/search?q=*",
            "https://www.aptanexus.com/api/search?query=",
            "https://www.aptanexus.com/api/aptamers?limit=10000",
            "https://aptamer-database.vercel.app/api/aptamers",
        ]
        for url in search_urls:
            raw_records = try_json_api(url)
            if raw_records:
                break

    if not raw_records:
        print("\n" + "=" * 60)
        print("Could not download data automatically.")
        print("")
        print("MANUAL FALLBACK:")
        print("1. Open https://www.aptanexus.com/ in your browser")
        print("2. Open DevTools (F12) → Network tab")
        print("3. Use the search bar or browse data")
        print("4. Look for API requests (XHR/Fetch) — copy the URL")
        print("5. Or try downloading directly:")
        print(f"   curl -o aptanexus.jsonl '{JSONL_URL}'")
        print("")
        print("If you get the file, place it as 'aptanexus_raw.jsonl'")
        print("in this directory and rerun with --parse-only flag")
        print("=" * 60)
        sys.exit(1)

    # Normalize records
    print(f"\nNormalizing {len(raw_records)} records...")
    normalized = [normalize_record(r) for r in raw_records]

    # Save
    print(f"\nSaving outputs...")
    save_outputs(normalized, raw_records)

    # Print summary
    print(f"\n{'=' * 60}")
    print(f"SUCCESS: {len(normalized)} aptamer records downloaded")
    print(f"{'=' * 60}")

    with_seq = sum(1 for r in normalized if r.get("sequence"))
    with_kd = sum(1 for r in normalized if r.get("kd_nM") is not None)
    with_target = sum(1 for r in normalized if r.get("target_name"))
    dna = sum(1 for r in normalized if r.get("aptamer_type") == "DNA")
    rna = sum(1 for r in normalized if r.get("aptamer_type") == "RNA")
    targets = len(set(r.get("target_name", "") for r in normalized if r.get("target_name")))
    tiers = {}
    for r in normalized:
        t = r.get("quality_tier", "unknown")
        tiers[t] = tiers.get(t, 0) + 1

    print(f"\nSummary:")
    print(f"  Total records:     {len(normalized)}")
    print(f"  With sequence:     {with_seq}")
    print(f"  With Kd value:     {with_kd}")
    print(f"  With target name:  {with_target}")
    print(f"  DNA aptamers:      {dna}")
    print(f"  RNA aptamers:      {rna}")
    print(f"  Unique targets:    {targets}")
    if tiers:
        print(f"  Quality tiers:     {tiers}")

    if with_kd > 0:
        kds = [r["kd_nM"] for r in normalized if r.get("kd_nM") is not None and r["kd_nM"] > 0]
        if kds:
            print(f"  Kd range:          {min(kds):.2f} – {max(kds):.2f} nM")

    # Print first record's keys for schema inspection
    if raw_records:
        print(f"\nRaw record schema (first record's keys):")
        first = raw_records[0]
        for k, v in sorted(first.items()):
            val_preview = str(v)[:60] if v is not None else "null"
            print(f"  {k}: {val_preview}")


if __name__ == "__main__":
    main()
