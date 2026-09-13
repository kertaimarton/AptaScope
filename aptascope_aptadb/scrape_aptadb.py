"""
AptaDB Scraper — Downloads all aptamer-target interaction records
from https://lmmd.ecust.edu.cn/aptadb/

Outputs:
  - aptadb_raw.json     (all records, structured)
  - aptadb_raw.csv      (flat table for quick inspection)

Usage:
  pip install requests beautifulsoup4 --break-system-packages
  python scrape_aptadb.py

Notes:
  - AptaDB has ~1350 interaction records (IDs 1 to ~1400, some gaps)
  - Script is rate-limited to 1 request/second to be respectful
  - Estimated runtime: ~25 minutes for full database
  - Saves progress every 50 records so you can resume if interrupted
"""

import requests
from bs4 import BeautifulSoup
import json
import csv
import time
import re
import os
import sys

BASE_URL = "https://lmmd.ecust.edu.cn/aptadb/result.php?id={}"
MAX_ID = 1400          # scan up to this ID (database has ~1350 entries)
DELAY = 1.0            # seconds between requests — don't hammer the server
CHECKPOINT_FILE = "aptadb_checkpoint.json"
OUTPUT_JSON = "aptadb_raw.json"
OUTPUT_CSV = "aptadb_raw.csv"
TIMEOUT = 15           # request timeout in seconds


def parse_result_page(html: str, interaction_id: int) -> dict | None:
    """
    Parse a single AptaDB result page into a structured record.
    Returns None if the page doesn't contain valid aptamer data.
    """
    soup = BeautifulSoup(html, "html.parser")

    # Check if page has actual content (some IDs are empty/404)
    tables = soup.find_all("table")
    if len(tables) < 2:
        return None

    record = {"interaction_id": interaction_id}

    # --- 1. Aptamer-target interaction table (first data table) ---
    # Look for the table with headers: Aptamer ID, Aptamer descriptor, etc.
    for table in tables:
        headers = [th.get_text(strip=True) for th in table.find_all("th")]
        rows = table.find_all("tr")

        # Interaction info table
        if "Aptamer ID" in headers and "Target name" in headers:
            for row in rows[1:]:  # skip header row
                cells = row.find_all("td")
                if len(cells) >= 7:
                    record["aptamer_id"] = cells[0].get_text(strip=True)
                    record["aptamer_descriptor"] = cells[1].get_text(strip=True)
                    record["target_chemistry"] = cells[2].get_text(strip=True)
                    record["target_name"] = cells[3].get_text(strip=True)
                    record["affinity_raw"] = cells[4].get_text(strip=True)
                    record["buffer_conditions"] = cells[5].get_text(strip=True)
                    # PubMed ID — extract from link or text
                    pubmed_link = cells[6].find("a")
                    if pubmed_link:
                        record["pubmed_id"] = pubmed_link.get_text(strip=True)
                        record["pubmed_url"] = pubmed_link.get("href", "")
                    else:
                        record["pubmed_id"] = cells[6].get_text(strip=True)

    # --- 2. Structure information ---
    # Look for the sequence and dot-bracket notation
    for table in tables:
        rows = table.find_all("tr")
        for row in rows:
            cells = row.find_all("td")
            if len(cells) == 2:
                label = cells[0].get_text(strip=True)
                value = cells[1].get_text(strip=True)

                if "Aptamer Sequence" in label and len(value) > 5:
                    record["sequence"] = value
                elif "optimal secondary structure" in label.lower():
                    record["mfe_structure"] = value
                elif "centroid secondary structure" in label.lower():
                    record["centroid_structure"] = value

    # --- 3. Aptamer detailed information ---
    # This table has a 4-column layout: Type | Detail | Type | Detail
    for table in tables:
        rows = table.find_all("tr")
        for row in rows:
            cells = row.find_all("td")
            # Process pairs of (label, value)
            cell_texts = [c.get_text(strip=True) for c in cells]

            for i in range(0, len(cell_texts) - 1, 2):
                label = cell_texts[i].lower()
                value = cell_texts[i + 1]

                if label == "aptamer chemistry":
                    record["aptamer_type"] = value  # DNA or RNA
                elif label == "length":
                    record["length_raw"] = value
                elif label == "gc content":
                    record["gc_content_raw"] = value
                elif label == "molecular weight":
                    record["molecular_weight"] = value
                elif "g-quadruplex" in label:
                    record["g_quadruplex_raw"] = value
                elif label == "g-score":
                    record["g_score"] = value
                elif label == "function":
                    record["function"] = value
                elif label == "applications":
                    record["applications"] = value

    # --- 4. Target information ---
    for table in tables:
        rows = table.find_all("tr")
        for row in rows:
            cells = row.find_all("td")
            cell_texts = [c.get_text(strip=True) for c in cells]

            for i in range(0, len(cell_texts) - 1, 2):
                label = cell_texts[i].lower()
                value = cell_texts[i + 1]

                if label == "uniprot id":
                    record["uniprot_id"] = value
                    # Also grab the link
                    links = cells[i + 1].find_all("a") if i + 1 < len(cells) else []
                    if links:
                        record["uniprot_url"] = links[0].get("href", "")
                elif label == "protein name":
                    record["protein_name"] = value
                elif label == "gene name(s)":
                    record["gene_names"] = value
                elif label == "organism":
                    record["organism"] = value
                elif label == "pdb id(s)":
                    record["pdb_ids"] = value

    # --- 5. Activity data ---
    # This section uses a different format (dl/list or inline)
    # Look for activity-related text
    text_content = soup.get_text()

    # Extract activity/affinity from the activity section
    activity_match = re.search(r"Activity\s+([\d.]+\s*[nμmuMpP]+)", text_content)
    if activity_match:
        record["activity_value"] = activity_match.group(1)

    # Extract assay description
    assay_section = soup.find(string=re.compile(r"Assay"))
    if assay_section:
        # Get the parent element and its next sibling text
        parent = assay_section.find_parent()
        if parent:
            next_text = parent.find_next_sibling()
            if next_text:
                record["assay_description"] = next_text.get_text(strip=True)[:500]

    # --- 6. Similar aptamers ---
    similar = []
    for table in tables:
        headers = [th.get_text(strip=True) for th in table.find_all("th")]
        if "Similarity" in headers:
            for row in table.find_all("tr")[1:]:
                cells = row.find_all("td")
                if len(cells) >= 4:
                    sim_record = {
                        "aptamer_id": cells[0].get_text(strip=True),
                        "chemistry": cells[1].get_text(strip=True),
                        "sequence": cells[2].get_text(strip=True),
                        "similarity": cells[3].get_text(strip=True),
                    }
                    # Get link to the similar aptamer
                    link = cells[0].find("a")
                    if link:
                        href = link.get("href", "")
                        id_match = re.search(r"id=(\d+)", href)
                        if id_match:
                            sim_record["interaction_id"] = int(id_match.group(1))
                    similar.append(sim_record)

    if similar:
        record["similar_aptamers"] = similar

    # --- 7. Target UniProt ID from the top of the page ---
    # "Target unique ID: O15243" pattern
    target_id_el = soup.find(string=re.compile(r"Target unique ID"))
    if target_id_el:
        parent = target_id_el.find_parent()
        if parent:
            link = parent.find("a")
            if link:
                record.setdefault("uniprot_id", link.get_text(strip=True))

    # --- Post-processing ---
    record = clean_record(record)

    # Validate: must have at least a sequence
    if "sequence" not in record or len(record.get("sequence", "")) < 5:
        return None

    return record


def clean_record(record: dict) -> dict:
    """Parse raw string values into clean types."""

    # Parse affinity to nanomolar
    affinity_raw = record.get("affinity_raw", "")
    if affinity_raw:
        record["kd_nM"] = parse_affinity_to_nM(affinity_raw)

    # Parse length to int
    length_raw = record.get("length_raw", "")
    length_match = re.search(r"(\d+)", length_raw)
    if length_match:
        record["length"] = int(length_match.group(1))

    # Parse GC content to float
    gc_raw = record.get("gc_content_raw", "")
    gc_match = re.search(r"([\d.]+)", gc_raw)
    if gc_match:
        record["gc_content"] = float(gc_match.group(1)) / 100.0  # convert % to fraction

    # Parse G-quadruplex
    gq_raw = record.get("g_quadruplex_raw", "")
    if "no" in gq_raw.lower() or "n/a" in gq_raw.lower():
        record["has_g_quadruplex"] = False
    elif gq_raw and "no" not in gq_raw.lower():
        record["has_g_quadruplex"] = True

    # Clean sequence — uppercase, strip whitespace
    if "sequence" in record:
        record["sequence"] = re.sub(r"[^ATGCUatgcu]", "", record["sequence"].upper())

    return record


def parse_affinity_to_nM(raw: str) -> float | None:
    """
    Convert affinity strings like '0.41 μM', '23 nM', '1.5 pM' to nanomolar.
    """
    raw = raw.strip()
    match = re.search(r"([\d.]+)\s*(nM|μM|uM|mM|pM|fM|M)", raw, re.IGNORECASE)
    if not match:
        return None

    value = float(match.group(1))
    unit = match.group(2).lower()

    conversions = {
        "fm": 1e-6,
        "pm": 1e-3,
        "nm": 1.0,
        "μm": 1e3,
        "um": 1e3,
        "mm": 1e6,
        "m": 1e9,
    }

    multiplier = conversions.get(unit)
    if multiplier is None:
        return None

    return round(value * multiplier, 4)


def load_checkpoint() -> tuple[list[dict], int]:
    """Load progress from checkpoint file."""
    if os.path.exists(CHECKPOINT_FILE):
        with open(CHECKPOINT_FILE, "r") as f:
            data = json.load(f)
            return data.get("records", []), data.get("last_id", 0)
    return [], 0


def save_checkpoint(records: list[dict], last_id: int):
    """Save progress to checkpoint file."""
    with open(CHECKPOINT_FILE, "w") as f:
        json.dump({"records": records, "last_id": last_id}, f)


def save_outputs(records: list[dict]):
    """Save final JSON and CSV outputs."""

    # JSON output (full data including similar aptamers)
    with open(OUTPUT_JSON, "w", encoding="utf-8") as f:
        json.dump(records, f, indent=2, ensure_ascii=False)
    print(f"  Saved {len(records)} records to {OUTPUT_JSON}")

    # CSV output (flat, without similar_aptamers nested field)
    if not records:
        return

    # Determine all flat keys (exclude nested fields)
    flat_keys = set()
    for r in records:
        for k, v in r.items():
            if not isinstance(v, (list, dict)):
                flat_keys.add(k)

    flat_keys = sorted(flat_keys)

    with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=flat_keys, extrasaction="ignore")
        writer.writeheader()
        for r in records:
            writer.writerow({k: r.get(k, "") for k in flat_keys})
    print(f"  Saved {len(records)} records to {OUTPUT_CSV}")


def scrape_all():
    """Main scraping loop."""
    records, last_id = load_checkpoint()

    if last_id > 0:
        print(f"Resuming from ID {last_id + 1} ({len(records)} records loaded)")
    else:
        print("Starting fresh scrape of AptaDB")

    start_id = last_id + 1
    session = requests.Session()
    session.headers.update({
        "User-Agent": "AptaScope-Research-Scraper/1.0 (academic research tool; contact: aptascope@example.com)",
        "Accept": "text/html",
    })

    empty_streak = 0  # track consecutive empty pages to detect end of database
    total_attempted = 0

    for interaction_id in range(start_id, MAX_ID + 1):
        total_attempted += 1

        try:
            url = BASE_URL.format(interaction_id)
            response = session.get(url, timeout=TIMEOUT)

            if response.status_code != 200:
                print(f"  ID {interaction_id}: HTTP {response.status_code}, skipping")
                empty_streak += 1
                if empty_streak > 30:
                    print(f"  30 consecutive misses — likely past end of database. Stopping.")
                    break
                time.sleep(DELAY)
                continue

            record = parse_result_page(response.text, interaction_id)

            if record:
                records.append(record)
                empty_streak = 0

                target = record.get("target_name", "unknown")
                kd = record.get("kd_nM", "?")
                seq_preview = record.get("sequence", "")[:30]
                print(f"  ID {interaction_id}: {record.get('aptamer_id', '?')} → {target} (Kd: {kd} nM) [{seq_preview}...]")
            else:
                empty_streak += 1
                if empty_streak > 30:
                    print(f"  30 consecutive empty pages — likely past end of database. Stopping.")
                    break

        except requests.exceptions.Timeout:
            print(f"  ID {interaction_id}: timeout, skipping")
            empty_streak += 1
        except requests.exceptions.ConnectionError:
            print(f"  ID {interaction_id}: connection error, waiting 10s and retrying...")
            time.sleep(10)
            try:
                response = session.get(BASE_URL.format(interaction_id), timeout=TIMEOUT)
                if response.status_code == 200:
                    record = parse_result_page(response.text, interaction_id)
                    if record:
                        records.append(record)
                        empty_streak = 0
            except Exception:
                print(f"  ID {interaction_id}: retry failed, skipping")
                empty_streak += 1
        except Exception as e:
            print(f"  ID {interaction_id}: unexpected error: {e}")
            empty_streak += 1

        # Checkpoint every 50 records
        if total_attempted % 50 == 0:
            save_checkpoint(records, interaction_id)
            print(f"  --- Checkpoint: {len(records)} records saved, at ID {interaction_id} ---")

        time.sleep(DELAY)

    # Final save
    save_checkpoint(records, interaction_id if 'interaction_id' in dir() else MAX_ID)
    save_outputs(records)

    print(f"\nDone. {len(records)} total aptamer records scraped.")
    print(f"Files: {OUTPUT_JSON}, {OUTPUT_CSV}")

    # Print summary stats
    if records:
        with_kd = sum(1 for r in records if r.get("kd_nM") is not None)
        dna = sum(1 for r in records if r.get("aptamer_type") == "DNA")
        rna = sum(1 for r in records if r.get("aptamer_type") == "RNA")
        targets = len(set(r.get("target_name", "") for r in records))
        print(f"\nSummary:")
        print(f"  Total records:  {len(records)}")
        print(f"  With Kd value:  {with_kd}")
        print(f"  DNA aptamers:   {dna}")
        print(f"  RNA aptamers:   {rna}")
        print(f"  Unique targets: {targets}")

        if with_kd > 0:
            kds = [r["kd_nM"] for r in records if r.get("kd_nM") is not None]
            print(f"  Kd range:       {min(kds):.2f} – {max(kds):.2f} nM")


if __name__ == "__main__":
    print("=" * 60)
    print("AptaDB Scraper for AptaScope")
    print("Target: https://lmmd.ecust.edu.cn/aptadb/")
    print(f"Scanning IDs 1–{MAX_ID}, delay {DELAY}s between requests")
    print("=" * 60)
    scrape_all()
