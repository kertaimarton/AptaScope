# AptaScope — Aptamer Sequence Analyzer & Design Explorer

## What this is

A single-page React web app that lets researchers analyze aptamer sequences against a database of ~2,500+ experimentally validated aptamer-target interactions. No backend — all computation runs in-browser. All data is embedded as a JSON module.

## Why it matters

Aptamers are short DNA/RNA sequences that bind specific targets (like antibodies, but synthetic and cheaper). Researchers design them through SELEX (a lab selection process), but there's no accessible tool to analyze a candidate sequence's biophysical features and compare it against known binders. Existing tools (AptaBERT, AptaBLE) are research code. Existing databases (AptaDB, UTexas) are browseable but not analyzable. AptaScope fills the gap.

---

## Build order

Complete each phase fully before moving to the next.

### Phase 1: Data acquisition and processing

**Goal:** Produce a single `aptamer_data.json` file with enriched aptamer records.

#### Step 1.1 — Download datasets

Download these two datasets:

1. **UTexas Aptamer Database** from Zenodo:
   - URL: `https://zenodo.org/records/8264921`
   - This is a ZIP containing CSV/Excel files with aptamer sequences, targets, Kd values, SELEX methods
   - Download and extract

2. **AptaDB bulk data** — if not directly downloadable as CSV, scrape the key fields from `https://lmmd.ecust.edu.cn/aptadb/` or use the supplementary data from their paper. The paper's supplementary materials (DOI: 10.1093/nar/gkae128) may have a downloadable table.

If either source is inaccessible (server down, paywall on download), proceed with whichever one works. The UTexas Zenodo dataset is the priority — it's the most reliably downloadable.

**Fallback:** If neither downloads cleanly, build a representative dataset by:
- Web-scraping Aptagen's Apta-Index (aptagen.com/apta-index) 
- Or using the Mendeley dataset at `https://data.mendeley.com/datasets/76jgjbgndr/1` (238 aptamers)
- Or hardcoding a curated dataset of ~200 well-known aptamers from literature (thrombin, VEGF, PDGF, cocaine, ATP aptamers are the classics — their sequences and Kd values are widely published)

The app must work with real data, not mock data.

#### Step 1.2 — Clean and unify

Python script (`process_data.py`) that:

1. Reads raw data from downloaded files
2. Extracts these fields per aptamer record:

```python
{
  "id": str,              # unique identifier
  "sequence": str,        # nucleotide sequence (uppercase, ATGCU only)
  "type": str,            # "DNA" or "RNA"  
  "length": int,          # len(sequence)
  "target_name": str,     # what it binds (e.g., "Thrombin", "VEGF")
  "target_type": str,     # "protein", "small_molecule", "cell", "other"
  "kd_nM": float | null,  # dissociation constant in nanomolar (lower = better binding)
  "selex_method": str | null,  # e.g., "CE-SELEX", "bead-based"
  "year": int | null,     # publication year
  "doi": str | null,      # paper DOI
  "source": str           # "utexas" or "aptadb"
}
```

3. Deduplicates on (sequence, target_name) — keep entry with more complete annotation
4. Normalizes Kd values to nanomolar (papers report in nM, μM, pM — convert all to nM)
5. Removes entries with no sequence

#### Step 1.3 — Compute biophysical features

Extend each record with computed features:

```python
{
  # ... base fields above, plus:
  "gc_content": float,       # (G+C) / length, range 0-1
  "a_freq": float,           # A count / length
  "t_freq": float,           # T (or U for RNA) count / length  
  "g_freq": float,           # G count / length
  "c_freq": float,           # C count / length
  "purine_ratio": float,     # (A+G) / length
  "complexity": float,       # Shannon entropy of dinucleotide frequencies (higher = more complex)
  "has_g_quadruplex": bool,  # regex match for G{3,}N{1,7}G{3,}N{1,7}G{3,}N{1,7}G{3,}
  "longest_repeat": int,     # longest homopolymer run (e.g., "AAAA" = 4)
  "predicted_mfe": float | null,  # minimum free energy from structure prediction (kcal/mol)
  "num_stems": int | null,   # predicted number of stem regions
  "num_loops": int | null    # predicted number of loop regions
}
```

**Shannon entropy calculation:**
```python
from collections import Counter
import math

def dinucleotide_entropy(seq):
    dinucs = [seq[i:i+2] for i in range(len(seq)-1)]
    counts = Counter(dinucs)
    total = sum(counts.values())
    probs = [c/total for c in counts.values()]
    return -sum(p * math.log2(p) for p in probs if p > 0)
```

**G-quadruplex detection:**
```python
import re
def has_g_quadruplex(seq):
    pattern = r'G{3,}.{1,7}G{3,}.{1,7}G{3,}.{1,7}G{3,}'
    return bool(re.search(pattern, seq))
```

**Secondary structure prediction:**
Try to install ViennaRNA (`pip install ViennaRNA`). If it installs:
```python
import RNA
def predict_structure(seq):
    structure, mfe = RNA.fold(seq)
    stems = structure.count('(')  # rough count
    loops = len(re.findall(r'\(\.+\)', structure))
    return mfe, stems, loops
```

If ViennaRNA doesn't install, implement a simplified Nussinov-style fold or skip the MFE/structure fields (set to null). The app should handle null values gracefully.

#### Step 1.4 — Output

Save as `src/data/aptamer_data.json` — an array of enriched records. Also save summary stats:

```python
{
  "total_records": int,
  "records_with_kd": int,
  "unique_targets": int,
  "dna_count": int,
  "rna_count": int,
  "kd_range": [min, max],
  "median_length": float,
  "median_gc": float,
  "generated_at": "ISO timestamp"
}
```

Save as `src/data/dataset_stats.json`.

---

### Phase 2: React application

**Stack:**
- React (functional components, hooks)
- Recharts for charts
- Tailwind CSS for styling
- No backend, no API calls, no localStorage

**Install and init:**
```bash
npx create-react-app aptascope --template default
cd aptascope
npm install recharts
# Tailwind setup via CDN in index.html is fine for this prototype
```

Copy `aptamer_data.json` and `dataset_stats.json` into `src/data/`.

#### App structure

```
src/
  data/
    aptamer_data.json
    dataset_stats.json
  components/
    Layout.jsx          # nav + main container
    Explorer.jsx        # searchable table view
    Analytics.jsx       # charts and stats dashboard
    Analyzer.jsx        # paste-a-sequence analyzer
    TargetLookup.jsx    # per-target deep dive
    SequenceDetail.jsx  # detail view for one aptamer
  utils/
    sequence.js         # sequence analysis functions (GC calc, entropy, etc.)
    similarity.js       # edit distance / alignment
    stats.js            # percentile ranking, distribution calcs
  App.jsx
  index.js
```

#### Navigation

Top bar with four tabs: **Explorer** | **Analytics** | **Analyzer** | **Targets**

#### View 1: Explorer

A searchable, filterable, sortable table of all aptamers.

**Filters (sidebar or top bar):**
- Text search (searches target name and sequence)
- Type toggle: All / DNA / RNA
- Kd range slider (log scale, from 0.01 nM to 100,000 nM)
- Length range slider (10–200 nt)
- Target type dropdown: All / Protein / Small Molecule / Cell / Other

**Table columns:**
- Target name
- Type (DNA/RNA badge)
- Sequence (truncated to 40 chars, monospace font, full on hover/click)
- Length
- Kd (nM) — color-coded: green <10nM, yellow 10-1000nM, red >1000nM
- GC%

**Click a row** → expand to show full detail: all computed features, link to DOI, sequence in copyable monospace block.

**Sort** by any column.

#### View 2: Analytics

Dashboard with 4-6 charts showing patterns across the whole dataset.

**Chart 1: GC Content vs. Binding Affinity**
- Scatter plot, x = GC content (0-1), y = log10(Kd) (inverted so better binding is higher)
- Color points by type (DNA = blue, RNA = orange)
- Add a trend line or moving average

**Chart 2: Length vs. Binding Affinity**
- Scatter plot, x = length, y = log10(Kd)
- Highlight the "sweet spot" zone if one emerges

**Chart 3: Affinity Distribution**
- Histogram of log10(Kd) values
- Show median and quartiles

**Chart 4: Top Targets**
- Horizontal bar chart of the 15 most-studied targets (by number of aptamers)
- Color bars by median Kd for that target

**Chart 5: Nucleotide Composition by Target Type**
- Grouped bar chart showing avg A/T/G/C frequencies for protein vs small-molecule vs cell targets

**Chart 6: Timeline**
- Line chart of aptamers published per year (if year data is available)

Each chart should have a one-line interpretation beneath it (e.g., "Most high-affinity binders cluster at 40-60% GC content").

#### View 3: Analyzer (core feature)

**Input section:**
- Large textarea: "Paste your aptamer sequence (DNA or RNA)"
- Auto-detect DNA vs RNA (presence of U = RNA, T = DNA)
- "Analyze" button
- Example buttons: pre-loaded sequences (thrombin aptamer, VEGF aptamer) so users can try it immediately

**Output section (appears after analysis):**

**Panel A: Computed Features**
A clean card grid showing:
- Length: X nt
- GC Content: X% (with percentile vs database: "72nd percentile")
- Sequence Complexity: X (Shannon entropy, with percentile)
- G-Quadruplex: Yes/No
- Longest Repeat: X nt
- Purine Ratio: X%
- MFE: X kcal/mol (if computed, with percentile)

Each metric should show a small inline sparkline or gauge showing where the input falls relative to the database distribution.

**Panel B: Radar Chart**
A radar/spider chart with axes: GC Content, Complexity, Length, Purine Ratio, (MFE if available). Plot two polygons:
1. The user's sequence (blue)
2. The median of high-affinity binders (Kd < 10nM) (green, semi-transparent)

This immediately shows how the input compares to known good aptamers.

**Panel C: Most Similar Known Aptamers**
Compute Levenshtein edit distance between the input sequence and every sequence in the database. Show top 5 matches:
- Table with: Target, Sequence (aligned), Kd, Similarity %, Source

**Panel D: Design Suggestions**
Based on the computed features, generate text suggestions:
- If GC < 0.35 or > 0.70: "GC content is outside the typical range for high-affinity binders (40-60%). Consider adjusting."
- If no structure predicted / very positive MFE: "Low structural stability. High-affinity aptamers typically form defined secondary structures."
- If G-quadruplex detected: "G-quadruplex motif detected. These structures are associated with strong target binding in many aptamer families."
- If length < 15 or > 80: "Sequence length is unusual. Most validated aptamers are 20-60 nt."

#### View 4: Target Lookup

**Input:** Dropdown/searchable select of all unique targets in the database.

**Output for selected target:**

- Number of known aptamers for this target
- Best Kd achieved
- Table of all aptamers for this target (sortable)
- "Design profile" — box showing the feature ranges of known binders for this target:
  - Optimal GC range
  - Optimal length range
  - Common structural motifs
  - Common nucleotide biases
- Comparison scatter: GC vs Kd for just this target's aptamers

---

### Phase 3: Design direction

This is a scientific tool, not a SaaS landing page. The design should feel like a well-made research instrument.

**Color palette:**
- Background: `#0D1117` (GitHub dark-style, not pure black)
- Surface/cards: `#161B22`
- Border: `#30363D`
- Text primary: `#E6EDF3`
- Text secondary: `#8B949E`
- Accent (primary actions, links): `#58A6FF` (clear blue)
- Accent success (high affinity): `#3FB950` (green)
- Accent warning (moderate): `#D29922` (amber)
- Accent danger (low affinity): `#F85149` (red)
- DNA color: `#58A6FF` (blue)
- RNA color: `#F0883E` (orange)

**Typography:**
- Use Inter for UI text (clean, readable, available from Google Fonts CDN)
- Use JetBrains Mono for sequences and code-like data (monospace, available from Google Fonts CDN)
- Body: 14px, headings scale: 18/24/32px
- Sequences always in monospace, with letter-spacing: 0.05em for readability

**Layout:**
- Max-width 1280px, centered
- Content left-aligned, not centered
- Cards with subtle borders, no heavy shadows
- Generous spacing between sections

**Key UI details:**
- Sequences should be displayed with nucleotide coloring: A=green, T=red, G=amber, C=blue (standard bioinformatics convention)
- Kd values should be formatted smartly: show as "4.2 nM" not "4.200000", use "< 1 nM" for sub-nanomolar
- Loading states for heavy computations (similarity search)
- Responsive but desktop-first (researchers use desktop)

---

### Phase 4: Finishing

After the app runs and all views work:

1. **Test with real sequences:**
   - Thrombin aptamer (HD1): `GGTTGGTGTGGTTGG` (Kd ~100 nM)
   - VEGF aptamer: `TGTGGGGGTGGACGGGCCGGGTAGA` (Kd ~1 nM)
   - ATP aptamer: `ACCTGGGGGAGTATTGCGGAGGAAGGT`

2. **Add an About section** (footer or modal):
   - What AptaScope is (one paragraph)
   - Data sources with links (UTexas, AptaDB)
   - What Kd means (one sentence)
   - Disclaimer: "This tool is for research exploration. Validate computationally-derived insights experimentally."

3. **CSV export button** in the Analyzer view — export computed features as CSV

4. **Copy sequence button** on every sequence display

---

## Technical constraints

- No backend, no database, no API calls to external services
- All data embedded in the app as JSON imports
- All computation in-browser (JavaScript)
- Levenshtein distance is O(n*m) per pair, and we run it against ~2,500 sequences — this may be slow. Optimize: pre-filter by length (skip sequences >2x or <0.5x the input length), and use a Web Worker if the UI freezes.
- Chart library: Recharts only (it's already available in React templates here)
- Style: Tailwind via CDN in index.html `<script>` tag, not PostCSS build
- Fonts: load Inter and JetBrains Mono from Google Fonts CDN in index.html

## File output

The final built app should be in `/home/claude/aptascope/` with:
- `npm start` runs the dev server
- All data files in `src/data/`
- All components working and connected
- Copy the build output or key files to `/mnt/user-data/outputs/` when complete
