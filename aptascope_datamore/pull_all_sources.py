#!/usr/bin/env python3
"""Pull and normalize new aptamer data sources for AptaScope.

Self-contained: outputs go to ./data_pulls/ only. Does NOT read or write
anything in the existing aptascope app repo, and does NOT dedup against
the existing 13,306-record dataset (that merge happens elsewhere, later).

Real endpoints below were confirmed live (2026-09-14) by inspecting each
site's actual JS bundle / network calls -- they differ from the guessed
URLs in the original spec doc for Ribocentre, AptCancerDB and AptBacterialDB.

Usage:
    python pull_all_sources.py                 # everything except ChEMBL (fast, minutes)
    python pull_all_sources.py --chembl-only    # ChEMBL enrichment only (slow, ~1-2 hours)
"""
import argparse
import json
import re
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import requests

HERE = Path(__file__).parent
OUT_DIR = HERE / "data_pulls"
OUT_DIR.mkdir(exist_ok=True)
DATAFIND_DIR = HERE / "aptascope_datafind"

session = requests.Session()
session.headers.update({"User-Agent": "Mozilla/5.0 (AptaScope research data collection)"})


def get_with_retry(url, attempts=4, backoff=2.0, **kwargs):
    kwargs.setdefault("timeout", 30)
    for attempt in range(1, attempts + 1):
        try:
            resp = session.get(url, **kwargs)
            resp.raise_for_status()
            return resp
        except requests.exceptions.RequestException as e:
            if attempt == attempts:
                raise
            wait = backoff * attempt
            print(f"  retry {attempt}/{attempts - 1} after error ({e}); waiting {wait:.0f}s")
            time.sleep(wait)

UNIT_TO_NM = {"pm": 1e-3, "nm": 1.0, "um": 1e3, "µm": 1e3, "μm": 1e3, "mm": 1e6, "m": 1e9}


def parse_affinity(raw) -> float | None:
    """Extract a Kd value in nM from a free-text affinity string."""
    if not raw:
        return None
    m = re.search(r"(\d+(?:\.\d+)?)\s*([pnuµμm]?M)\b", str(raw), re.IGNORECASE)
    if not m:
        return None
    value, unit = float(m.group(1)), m.group(2).lower()
    factor = UNIT_TO_NM.get(unit)
    if factor is None:
        return None
    nm = value * factor
    return nm if 1e-6 < nm < 1e12 else None


def detect_type(raw) -> str | None:
    if not raw:
        return None
    s = str(raw).lower()
    if "rna" in s:
        return "RNA"
    if "dna" in s:
        return "DNA"
    return None


def to_int(raw) -> int | None:
    try:
        return int(str(raw).strip())
    except (ValueError, TypeError):
        return None


def clean_seq(raw) -> str | None:
    if not raw:
        return None
    seq = re.sub(r"[^ATGCUatgcu]", "", str(raw).upper())
    return seq if len(seq) >= 5 else None


# --------------------------------------------------------------------------
# Source pullers (raw)
# --------------------------------------------------------------------------

def pull_ribocentre():
    """Static JSON data file backing the /sequences/ search page (not a REST API)."""
    url = "https://aptamer.ribocentre.org/apidata/sequences_cleaned.json"
    resp = get_with_retry(url, timeout=30)
    return resp.json()["Sheet1"]


def pull_utexas():
    """Zenodo record 8264921 concept redirects to the latest version; grab the newer xlsx."""
    api = "https://zenodo.org/api/records/8264921"
    resp = get_with_retry(api, timeout=30, allow_redirects=True)
    files = resp.json().get("files", [])
    target = next((f for f in files if "Sept2023" in f.get("key", "")), files[0])
    file_resp = get_with_retry(target["links"]["self"], timeout=120)
    xlsx_path = OUT_DIR / "utexas_source.xlsx"
    xlsx_path.write_bytes(file_resp.content)
    df = pd.read_excel(xlsx_path, sheet_name="Dataset")
    return df.to_dict(orient="records")


def pull_aptcancerdb():
    """React SPA; real data lives at /api/aptamers (server ignores the pageSize cap
    and returns everything, but the `total` field is authoritative -- loop pages
    defensively in case that changes)."""
    records = []
    page = 1
    while True:
        url = f"https://webs.iiitd.edu.in/raghava/aptcancerdb/api/aptamers?page={page}&pageSize=1500"
        resp = get_with_retry(url, timeout=90)
        data = resp.json()
        rows = data.get("rows", [])
        if not rows:
            break
        records.extend(rows)
        total = data.get("total", len(records))
        print(f"  aptcancerdb: {len(records)}/{total}")
        if len(records) >= total:
            break
        page += 1
        time.sleep(0.3)
    return records


def pull_aptbacterialdb():
    """PHP backend at /api/query.php, hard-capped at 100 records per page.
    Writes a partial checkpoint after every page so a mid-run connection drop
    (this host has been flaky) doesn't lose already-fetched pages on retry."""
    checkpoint_path = OUT_DIR / "aptbacterialdb_checkpoint.json"
    records = []
    offset, limit = 0, 100
    if checkpoint_path.exists():
        records = json.loads(checkpoint_path.read_text())
        offset = len(records)
        print(f"  resuming from checkpoint at offset={offset}")
    while True:
        url = f"https://webs.iiitd.edu.in/raghava/aptbacterialdb/api/query.php?limit={limit}&offset={offset}"
        resp = get_with_retry(url, timeout=30)
        rows = resp.json().get("data", [])
        if not rows:
            break
        records.extend(rows)
        checkpoint_path.write_text(json.dumps(records))
        print(f"  aptbacterialdb: {len(records)} (offset={offset})")
        offset += limit
        time.sleep(0.3)
    checkpoint_path.unlink(missing_ok=True)
    return records


def load_local_aptagen():
    return json.loads((DATAFIND_DIR / "aptagen_raw.json").read_text())


def load_local_rnaapt3d():
    return json.loads((DATAFIND_DIR / "rnaapt3d_raw.json").read_text())


# --------------------------------------------------------------------------
# Normalizers -> common schema
# (mirrors the field names already used in the app's aptamer_data.json so
# the later merge is a straight concat, not a remapping exercise)
# --------------------------------------------------------------------------

def normalize_ribocentre(records):
    out = []
    for r in records:
        seq = clean_seq(r.get("Sequence"))
        if not seq:
            continue
        out.append({
            "id": f"RIBO-{r.get('ID', '')}",
            "sequence": seq,
            "type": "RNA" if "U" in seq else "DNA",
            "length": r.get("Length") or len(seq),
            "target_name": (r.get("Ligand") or "").strip() or None,
            "target_type": r.get("Type"),
            "kd_nM": parse_affinity(r.get("Affinity")),
            "affinity_raw": r.get("Affinity"),
            "year": r.get("Year"),
            "doi": r.get("Link to PubMed Entry"),
            "aptamer_name": r.get("Named") or r.get("Article name"),
            "gc_content": r.get("GC Content"),
            "ligand_description": r.get("Ligand Description"),
            "ligand_id": r.get("Ligand Information \n(CAS/sequence)"),
            "source": "ribocentre",
        })
    return out


def normalize_utexas(records):
    out = []
    for i, r in enumerate(records):
        seq = clean_seq(r.get("Aptamer Sequence"))
        if not seq:
            continue
        year = r.get("Year of Paper")
        out.append({
            "id": f"UTEX-{i+1:05d}",
            "sequence": seq,
            "type": detect_type(r.get("Type of Nucleic Acid")) or ("RNA" if "U" in seq else "DNA"),
            "length": r.get("Sequence Length") or len(seq),
            "target_name": (str(r.get("Target ", "")).strip() or None),
            "kd_nM": r.get("Kd (nM)") if pd.notna(r.get("Kd (nM)")) else parse_affinity(r.get("Affinity")),
            "affinity_raw": r.get("Affinity"),
            "year": int(year) if pd.notna(year) else None,
            "doi": r.get("Journal DOI"),
            "aptamer_name": r.get("Name of Aptamer"),
            "gc_content": r.get("GC Content "),
            "selex_method": r.get("Pool Type"),
            "citation": r.get("Citation"),
            "source": "utexas",
        })
    return out


def normalize_aptcancerdb(records):
    out = []
    for r in records:
        seq = clean_seq(r.get("sequence"))
        if not seq:
            continue
        out.append({
            "id": f"ACDB-{r.get('id', '')}",
            "sequence": seq,
            "type": "RNA" if "U" in seq else "DNA",
            "length": to_int(r.get("length")) or len(seq),
            "target_name": (r.get("target") or "").strip() or None,
            "target_type": "protein",
            "kd_nM": parse_affinity(r.get("kdValue")),
            "affinity_raw": r.get("kdValue"),
            "year": to_int(r.get("year")),
            "doi": r.get("pmid"),
            "aptamer_name": r.get("name"),
            "cancer_type": r.get("cancerType"),
            "cell_line": r.get("cellLine"),
            "application": r.get("activityRole"),
            "modification": r.get("modification"),
            "selex_method": r.get("method"),
            "objective": r.get("objectiveMechanism"),
            "source": "aptcancerdb",
        })
    return out


def normalize_aptbacterialdb(records):
    out = []
    for r in records:
        seq = clean_seq(r.get("Sequence (5′ to 3′)"))
        if not seq:
            continue
        out.append({
            "id": f"ABDB-{r.get('AptBacDB_ID', '')}",
            "sequence": seq,
            "type": detect_type(r.get("Type")) or ("RNA" if "U" in seq else "DNA"),
            "length": to_int(r.get("Length")) or len(seq),
            "target_name": (r.get("Target") or "").strip() or None,
            "target_type": "protein",
            "kd_nM": parse_affinity(r.get("Kd  value")),
            "affinity_raw": r.get("Kd  value"),
            "year": to_int(r.get("Year")),
            "doi": r.get("PMID/DOI"),
            "aptamer_name": r.get("Aptamer name"),
            "bacterial_species": r.get("Target Organism"),
            "application": r.get("Activity Role"),
            "modification": r.get("Modification"),
            "selex_method": r.get("Method"),
            "objective": r.get("Objective/Mechanism"),
            "source": "aptbacterialdb",
        })
    return out


def normalize_aptagen(records):
    out = []
    for i, r in enumerate(records):
        seq = clean_seq(r.get("sequence"))
        if not seq:
            continue  # ~8% of Aptagen listings are commercial entries with no disclosed sequence
        out.append({
            "id": f"APTAGEN-{i+1:05d}",
            "sequence": seq,
            "type": r.get("aptamer_type") or ("RNA" if "U" in seq else "DNA"),
            "length": len(seq),
            "target_name": r.get("target_name"),
            "target_type": r.get("target_type"),
            "kd_nM": r.get("kd_nM"),
            "affinity_raw": r.get("affinity_raw"),
            "aptamer_name": r.get("name"),
            "aptagen_id": r.get("aptagen_id"),
            "source": "aptagen",
        })
    return out


def normalize_rnaapt3d(records):
    out = []
    for r in records:
        seq = clean_seq(r.get("sequence"))
        if not seq:
            continue
        target_raw = r.get("target") or ""
        out.append({
            "id": f"RNAAPT3D-{r.get('monad_id', '')}",
            "sequence": seq,
            "type": r.get("aptamer_type") or "RNA",
            "length": r.get("length") or len(seq),
            "target_name": re.sub(r"\s*\([^)]*\)\s*$", "", target_raw).strip() or None,
            "uniprot_id": r.get("uniprot_id"),
            "structure_pdb_urls": [s.get("pdb_url") for s in r.get("structures", []) if s.get("pdb_url")],
            "dot_bracket": (r.get("structures") or [{}])[0].get("dot_bracket"),
            "source": "rnaapt3d",
        })
    return out


# --------------------------------------------------------------------------
# Enrichment
# --------------------------------------------------------------------------

def enrich_uniprot(uniprot_ids):
    enriched = {}
    ids = sorted(set(uniprot_ids))
    for i, uid in enumerate(ids):
        try:
            resp = session.get(f"https://rest.uniprot.org/uniprotkb/{uid}.json", timeout=15,
                                headers={"Accept": "application/json"})
            if resp.status_code != 200:
                continue
            data = resp.json()
            enriched[uid] = {
                "uniprot_id": uid,
                "protein_name": data.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", ""),
                "gene_name": next((g.get("geneName", {}).get("value", "") for g in data.get("genes", [])), ""),
                "organism": data.get("organism", {}).get("scientificName", ""),
                "function": next((c.get("texts", [{}])[0].get("value", "") for c in data.get("comments", []) if c.get("commentType") == "FUNCTION"), ""),
                "subcellular_location": [loc.get("location", {}).get("value", "") for c in data.get("comments", []) if c.get("commentType") == "SUBCELLULAR LOCATION" for loc in c.get("subcellularLocations", [])],
                "diseases": [d.get("disease", {}).get("diseaseId", "") for c in data.get("comments", []) if c.get("commentType") == "DISEASE" for d in [c]],
                "molecular_weight": data.get("sequence", {}).get("molWeight"),
                "length_aa": data.get("sequence", {}).get("length"),
            }
        except Exception as e:
            print(f"  {uid}: {e}")
        if (i + 1) % 10 == 0:
            print(f"  UniProt: {i+1}/{len(ids)}")
        time.sleep(0.3)
    return enriched


def enrich_pdb_aptamer_structures():
    """No source pulled here carries real RCSB PDB codes (RNAapt3D hosts its own
    predicted structures, not RCSB entries). Instead of a per-ID lookup that would
    enrich nothing, pull RCSB's full-text search for all aptamer-related deposited
    structures as a standalone reference table."""
    search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
    query = {
        "query": {"type": "terminal", "service": "full_text", "parameters": {"value": "aptamer"}},
        "return_type": "entry",
        "request_options": {"results_content_type": ["experimental"], "return_all_hits": True},
    }
    resp = session.post(search_url, json=query, timeout=30)
    resp.raise_for_status()
    pdb_ids = [hit["identifier"] for hit in resp.json().get("result_set", [])]
    print(f"  Found {len(pdb_ids)} aptamer-related PDB entries")

    enriched = {}
    for i, pdb_id in enumerate(pdb_ids):
        try:
            resp = session.get(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}", timeout=15)
            if resp.status_code != 200:
                continue
            data = resp.json()
            enriched[pdb_id] = {
                "pdb_id": pdb_id,
                "title": data.get("struct", {}).get("title", ""),
                "method": (data.get("exptl") or [{}])[0].get("method", ""),
                "resolution": (data.get("rcsb_entry_info", {}).get("resolution_combined") or [None])[0],
                "deposit_date": data.get("rcsb_accession_info", {}).get("deposit_date", ""),
                "polymer_count": data.get("rcsb_entry_info", {}).get("polymer_entity_count", 0),
                "viewer_url": f"https://www.rcsb.org/3d-view/{pdb_id}",
            }
        except Exception:
            pass
        if (i + 1) % 50 == 0:
            print(f"  PDB: {i+1}/{len(pdb_ids)}")
        time.sleep(0.2)
    return enriched


_RETRY = object()  # sentinel: transient failure, don't checkpoint it -- retry next run


def _lookup_chembl_one(target):
    """One target -> (target, result). result is a dict on a real match,
    False on a confirmed no-match (checkpointed either way so it's never
    re-queried again), or _RETRY on a transient error (left out of the
    checkpoint so the next run retries it). Uses its own Session since
    requests.Session isn't safe to share across threads."""
    base = "https://www.ebi.ac.uk/chembl/api/data"
    local_session = requests.Session()
    local_session.headers.update(session.headers)
    try:
        resp = local_session.get(f"{base}/target/search.json", params={"q": target, "limit": 1}, timeout=15)
        if resp.status_code != 200:
            return target, _RETRY
        targets = resp.json().get("targets", [])
        if not targets:
            return target, False

        chembl_target = targets[0]
        chembl_id = chembl_target.get("target_chembl_id", "")
        drug_resp = local_session.get(f"{base}/mechanism.json",
                                       params={"target_chembl_id": chembl_id, "limit": 10}, timeout=15)
        drugs = []
        if drug_resp.status_code == 200:
            drugs = [
                {"drug_name": m.get("molecule_chembl_id", ""),
                 "mechanism": m.get("mechanism_of_action", ""),
                 "action_type": m.get("action_type", "")}
                for m in drug_resp.json().get("mechanisms", [])
            ]
        return target, {
            "chembl_id": chembl_id,
            "target_type": chembl_target.get("target_type", ""),
            "organism": chembl_target.get("organism", ""),
            "drugs": drugs,
        }
    except Exception:
        return target, _RETRY


def enrich_chembl(target_names, workers=12):
    checkpoint_path = OUT_DIR / "chembl_checkpoint.json"
    enriched = {}
    if checkpoint_path.exists():
        raw = json.loads(checkpoint_path.read_text())
        enriched = {k: (v if v else False) for k, v in raw.items()}
        print(f"  resuming from checkpoint: {len(enriched)} targets already attempted")

    names = sorted({n.strip() for n in target_names if n and n.strip()})
    todo = [n for n in names if n not in enriched]
    print(f"  {len(names)} unique targets total, {len(todo)} left to attempt", flush=True)

    completed = 0
    lock = threading.Lock()
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(_lookup_chembl_one, t) for t in todo]
        for future in as_completed(futures):
            target, result = future.result()
            if result is _RETRY:
                continue
            with lock:
                enriched[target] = result
                completed += 1
                if completed % 100 == 0:
                    checkpoint_path.write_text(json.dumps(enriched))
                    print(f"  ChEMBL: {completed}/{len(todo)}", flush=True)

    checkpoint_path.write_text(json.dumps(enriched))
    hits = {k: v for k, v in enriched.items() if v}
    print(f"  Done: {len(hits)}/{len(enriched)} targets matched in ChEMBL "
          f"({len(names) - len(enriched)} left needing a retry pass)")
    if len(enriched) >= len(names):
        checkpoint_path.unlink(missing_ok=True)
    return hits


# --------------------------------------------------------------------------
# Combine
# --------------------------------------------------------------------------

def combine_and_dedup(all_normalized):
    seen = set()
    combined = []
    dupes = 0
    for r in all_normalized:
        key = (r["sequence"], (r.get("target_name") or "").strip().lower())
        if key in seen:
            dupes += 1
            continue
        seen.add(key)
        # Real SELEX aptamers are essentially never this short; these are most likely
        # truncated/mismapped fields from source scrapes (concentrated in Aptagen).
        # Flagged rather than dropped -- decide at merge time, don't destroy data here.
        if len(r["sequence"]) < 15:
            r["quality_flag"] = "short_sequence"
        combined.append(r)
    flagged = sum(1 for r in combined if r.get("quality_flag") == "short_sequence")
    print(f"Combined: {len(combined)} unique records ({dupes} duplicates removed, "
          f"dedup key = sequence + target_name, within these new sources only)")
    print(f"  {flagged} records flagged quality_flag=short_sequence (<15nt) -- kept, not dropped")
    return combined


def save(name, data):
    path = OUT_DIR / name
    with open(path, "w") as f:
        json.dump(data, f, indent=2, ensure_ascii=False, default=str)
    print(f"  -> {path} ({len(data) if hasattr(data, '__len__') else '?'} entries)")


# --------------------------------------------------------------------------

def cached_pull(filename, pull_fn, label):
    path = OUT_DIR / filename
    if path.exists():
        print(f"=== {label} (cached) ===")
        return json.loads(path.read_text())
    print(f"=== {label} ===")
    raw = pull_fn()
    save(filename, raw)
    return raw


def run_main_pull():
    ribo_raw = cached_pull("ribocentre_raw.json", pull_ribocentre, "SOURCE 1: Ribocentre")
    ribo = normalize_ribocentre(ribo_raw)

    utexas_raw = cached_pull("utexas_raw.json", pull_utexas, "SOURCE 2: UTexas")
    utexas = normalize_utexas(utexas_raw)

    acdb_raw = cached_pull("aptcancerdb_raw.json", pull_aptcancerdb, "SOURCE 3: AptCancerDB")
    acdb = normalize_aptcancerdb(acdb_raw)

    abdb_raw = cached_pull("aptbacterialdb_raw.json", pull_aptbacterialdb, "SOURCE 4: AptBacterialDB")
    abdb = normalize_aptbacterialdb(abdb_raw)

    print("=== SOURCE 5: Aptagen (local, already pulled) ===")
    aptagen = normalize_aptagen(load_local_aptagen())

    print("=== SOURCE 6: RNAapt3D (local, already pulled) ===")
    rnaapt3d = normalize_rnaapt3d(load_local_rnaapt3d())

    for name, recs in [("ribocentre", ribo), ("utexas", utexas), ("aptcancerdb", acdb),
                        ("aptbacterialdb", abdb), ("aptagen", aptagen), ("rnaapt3d", rnaapt3d)]:
        print(f"  normalized {name}: {len(recs)} records")

    print("=== COMBINE ===")
    combined = combine_and_dedup(ribo + utexas + acdb + abdb + aptagen + rnaapt3d)
    save("aptascope_expansion.json", combined)

    print("=== ENRICHMENT: UniProt ===")
    uniprot_ids = {r["uniprot_id"] for r in combined if r.get("uniprot_id")}
    uniprot_data = enrich_uniprot(uniprot_ids)
    save("uniprot_enrichment.json", uniprot_data)

    print("=== ENRICHMENT: PDB (aptamer structures reference table) ===")
    pdb_data = enrich_pdb_aptamer_structures()
    save("pdb_enrichment.json", pdb_data)

    stats = {
        "generated": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "source_counts": {
            "ribocentre": len(ribo), "utexas": len(utexas), "aptcancerdb": len(acdb),
            "aptbacterialdb": len(abdb), "aptagen": len(aptagen), "rnaapt3d": len(rnaapt3d),
        },
        "combined_total": len(combined),
        "uniprot_enriched": len(uniprot_data),
        "pdb_structures_found": len(pdb_data),
        "chembl_enriched": "run with --chembl-only separately",
    }
    save("pull_summary.json", stats)
    print(json.dumps(stats, indent=2))
    print("\nDone. Run `python pull_all_sources.py --chembl-only` next (slow, background-friendly).")


def run_chembl_only():
    combined_path = OUT_DIR / "aptascope_expansion.json"
    if not combined_path.exists():
        print("Run the main pull first (no --chembl-only flag) -- aptascope_expansion.json not found.")
        sys.exit(1)
    combined = json.loads(combined_path.read_text())
    target_names = {r.get("target_name") for r in combined if r.get("target_name")}
    print(f"=== ENRICHMENT: ChEMBL ({len(target_names)} unique targets) ===")
    chembl_data = enrich_chembl(target_names)
    save("chembl_enrichment.json", chembl_data)
    print(f"Done. {len(chembl_data)}/{len(target_names)} targets matched in ChEMBL.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--chembl-only", action="store_true")
    args = parser.parse_args()
    if args.chembl_only:
        run_chembl_only()
    else:
        run_main_pull()
