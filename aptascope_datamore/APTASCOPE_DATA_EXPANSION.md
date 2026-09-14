# AptaScope Data Expansion — Merge Instructions

## Status

Data collection is **done**. Merging it into the main AptaScope app is **not done yet** — that's what this doc is for.

This folder (`aptascope_datamore`) was built and run with no access to the main AptaScope app repo (the one containing `src/data/aptascope_merged.json`, currently 13,306 records from AptaNexus + AptaDB). Everything here was pulled and normalized independently, deduplicated only against itself. If you're reading this after copying this folder somewhere new — find the main app repo (search for `aptascope_merged.json`) and do the merge described below.

## What's in `data_pulls/`

| File | What it is |
|---|---|
| `aptascope_expansion.json` | **The main deliverable.** 5,758 normalized aptamer records from 6 sources, deduped against each other only |
| `uniprot_enrichment.json` | 25 UniProt target records (function, organism, disease), keyed by UniProt ID |
| `pdb_enrichment.json` | 437 aptamer-related RCSB PDB structures (title, method, resolution), keyed by PDB ID |
| `chembl_enrichment.json` | 1,779 targets with known drug/mechanism data, keyed by target name |
| `pull_summary.json` | run stats + the data-quality caveats below |
| `*_raw.json`, `utexas_source.xlsx` | untouched raw pulls per source, kept for provenance |

`pull_all_sources.py` (repo root, next to this file) is the script that produced everything above. It's checkpointed and safe to re-run — cached raw pulls are reused, and `--chembl-only` resumes the ChEMBL pass without re-querying already-resolved targets.

### Source breakdown (post-dedup, in `aptascope_expansion.json`)

| Source | Records | Real endpoint used (differs from any older notes guessing at URLs) |
|---|---|---|
| AptBacterialDB | 1,578 | `webs.iiitd.edu.in/raghava/aptbacterialdb/api/query.php` — real paginated JSON API |
| UTexas | 1,344 | Zenodo record 8264921 → xlsx download |
| AptCancerDB | 1,316 | `webs.iiitd.edu.in/raghava/aptcancerdb/api/aptamers` — real React SPA backend API |
| Aptagen | 964 | local copy from `aptascope_datafind/aptagen_raw.json`; ~8% of listings had no disclosed sequence and were excluded |
| Ribocentre | 507 | `aptamer.ribocentre.org/apidata/sequences_cleaned.json` — a static JSON file, not a REST API |
| RNAapt3D | 49 | local copy from `aptascope_datafind/rnaapt3d_raw.json`; the only source here with real UniProt IDs attached |
| **Total** | **5,758** | 1,089 cross-source duplicates already removed (dedup key: sequence + target_name) |

### Record schema

Every record has: `id, sequence, type, length, target_name, kd_nM, source`. Fields beyond that vary by source (`cancer_type` on aptcancerdb, `bacterial_species` on aptbacterialdb, `uniprot_id` on rnaapt3d, etc.) — inspect a sample record per source or check `pull_all_sources.py`'s `normalize_*` functions.

**Notably absent:** the derived structural/compositional features the main app computes for every existing record — `gc_content` (present only on ribocentre/utexas here), `a_freq`, `t_freq`, `g_freq`, `c_freq`, `purine_ratio`, `complexity`, `has_g_quadruplex`, `longest_repeat`, `predicted_mfe`, `num_stems`, `num_loops`. See step 2 below.

### Known data-quality flags (deliberately left in, not fixed — your call)

- **107 records** carry `quality_flag: "short_sequence"` (<15nt), concentrated in Aptagen. Real SELEX aptamers are essentially never this short, so these are most likely truncated/mismapped fields from that scrape — but could include a legitimate minimal-motif aptamer here or there. Kept rather than dropped; decide at merge time.
- **24 target names** never resolved in ChEMBL (`A1`, `H2`, `71`, `PB`, etc. — see `pull_summary.json`). These are leftover well-plate/clone codes from the same Aptagen scrape, not real target names. Not worth re-attempting.

## Before merging into the main app

1. **Check for overlap on UTexas and Ribocentre specifically.** As of this pull, the main app repo's `src/data/aptamer_data.json` already contained 1,757 normalized records (1,415 UTexas + 342 Ribocentre) pulled independently in an earlier session, sitting unmerged into `aptascope_merged.json`. This directory's `aptascope_expansion.json` *also* contains fresh UTexas (1,344) and Ribocentre (507) pulls of the same two source databases, done separately with different normalization — counts differ because of different filtering, not different underlying data. **Don't add both sets** — reconcile them (dedupe one against the other on sequence, or just pick the more complete/feature-rich version) before folding into the merged dataset. AptCancerDB, AptBacterialDB, Aptagen, and RNAapt3D have no such prior overlap — those are safe to add wholesale (after the global dedup pass in step 3).

2. **Recompute structural features for every new record before merging.** The main app's `merge_databases.py` has `compute_features()` / `enrich()` functions that compute `gc_content`, `complexity`, `has_g_quadruplex`, `predicted_mfe`, `num_stems`, `num_loops`, etc. for every record (uses ViennaRNA for structure prediction). Reuse those functions rather than reimplementing them — run every sequence from `aptascope_expansion.json` through the app's existing pipeline so the merged dataset has a consistent schema. Explorer, Analytics, and `dataset_stats.json` all depend on these fields being present on every record.

3. **Dedupe against the full existing dataset, not just internally.** `aptascope_expansion.json` is only deduped against itself. Use the same key convention as `merge_databases.py`'s `deduplicate()` (sequence + canonicalized target name) against both `aptascope_merged.json` (13,306) and `aptamer_data.json` (1,757).

4. **Load the three enrichment files as lookup tables, not aptamer records** — they're keyed by UniProt ID / PDB ID / target name, meant to be joined at render time, not concatenated into the aptamer array.

5. **Regenerate `dataset_stats.json`** after merging — same shape as the existing file (`total_records`, `source_counts`, `kd_range`, quality tiers, etc.).

## App UI work needed after the data merge

- **`TargetLookup.jsx`** already has a `resolveUniProtId()` helper and a conditional UniProt block (`{uniprotId && (...)}`) — extend it with `function` / `organism` / `diseases` from `uniprot_enrichment.json`, and add a "Known drugs for this target" section sourced from `chembl_enrichment.json`.
- **`Explorer.jsx`** — add a "Has 3D Structure" filter badge (cross-reference `pdb_enrichment.json` by target name; none of these 6 sources carry a real per-record RCSB PDB ID to join on directly), plus "Cancer Type" (aptcancerdb) and "Bacterial Species" (aptbacterialdb) columns.
- **About / Data Sources section** — update source counts and credit the 4 new databases + 3 enrichment APIs.

## Expected numbers after a full merge

| Stage | Records |
|---|---|
| Current live app (`aptascope_merged.json`) | 13,306 |
| Already-pulled, unmerged UTexas + Ribocentre (`aptamer_data.json`) | 1,757 (reconcile against the row below first — same two sources, don't double-count) |
| This directory's new sources, all 6 combined | 5,758 (of which AptCancerDB + AptBacterialDB + Aptagen + RNAapt3D = 3,907 are unambiguously new) |
| **Realistic final total after cross-dedup** | **roughly 17,500–18,500** — exact number depends on how much overlap step 1 finds between the two UTexas/Ribocentre pulls |
