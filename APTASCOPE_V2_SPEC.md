# AptaScope v2 — Feature Expansion Spec

## Current state

AptaScope v1 is running with a React app containing four views: Explorer, Analytics, Analyzer, and Target Lookup. It currently uses a partial dataset. This spec covers:

1. Merging two complete databases into the app
2. Six new feature branches, ordered by cost (free builds first)

---

## PHASE 0: Database Merge

### What we have

**AptaNexus** (primary — largest open aptamer dataset that exists):
- File: `aptanexus_raw.jsonl` — 12,535 records
- 100% have sequences, 99.97% have target names
- 51.5% have binding affinity values
- 12,531 have secondary structure (dot-bracket) and MFE
- 1,875 unique targets
- Quality tiers: A (3,421), B (846), C (6,874), P (1,394)
- Field mapping:

```
"Aptamer sequence"    → sequence
"Target name"         → target_name
"Target type"         → target_type  (e.g. "1. Proteins & Peptides")
"Affinity"            → affinity_raw (string — needs parsing)
"MFE"                 → predicted_mfe (float, kcal/mol, already numeric)
"SecStr_dotbracket"   → mfe_structure (dot-bracket notation)
"Level"               → quality_tier  (A/B/C/P)
"Year"                → year
"Doi"                 → doi
"Journal"             → journal
"Article title"       → article_title
"External_ID"         → external_id (usually UniProt ID)
"External_Name"       → external_name
"Gene_Symbol"         → gene_symbol
"ID_Type"             → id_type
"pKd"                 → pkd (float, -log10 of Kd in M — already computed)
"Buffer condition"    → buffer
"Best"                → is_best (boolean — best aptamer for that target)
"Sequence ID"         → sequence_id (name given in original paper)
"SecStr_image"        → structure_image_path (SVG path, not useful offline)
```

**Affinity parsing notes:** The `Affinity` field is messy. Examples:
- `"24 nM"` → simple
- `"24 nM (SPR),12 nM (MST),1.81 nM (radiolabelling filter assay)"` → multiple methods, take the best (lowest)
- `"N/A"` → null
- `"0.05-0.1 μM"` → range, take midpoint and convert

When `pKd` is present and `Affinity` is hard to parse, derive Kd from pKd:
```
Kd (nM) = 10^(-pKd) * 1e9
```

**Target type cleanup:** Strip the leading number prefix:
- `"1. Proteins & Peptides"` → `"Protein"`
- `"7. Small Molecules & Ions"` → `"Small Molecule"`
- `"3. Cells, Tissues & Subcellular Structures"` → `"Cell"`
- `"2. Nucleic Acids"` → `"Nucleic Acid"`
- `"4. Microorganisms & Viruses"` → `"Microorganism"`
- Everything else → `"Other"`

**AptaDB** (secondary — richer per-record metadata):
- File: `aptadb_raw.json` — ~1,350 records (from scraper)
- Has detailed assay descriptions, buffer conditions, similar aptamers
- Field mapping (already matches the v1 schema from scrape_aptadb.py):

```
"aptamer_id"          → aptamer_id
"sequence"            → sequence
"target_name"         → target_name
"target_chemistry"    → target_type
"affinity_raw"        → affinity_raw
"kd_nM"               → kd_nM (already parsed to nM)
"aptamer_type"        → aptamer_type (DNA/RNA)
"length"              → length
"gc_content"          → gc_content (fraction 0-1)
"mfe_structure"       → mfe_structure (dot-bracket)
"centroid_structure"  → centroid_structure
"molecular_weight"    → molecular_weight
"has_g_quadruplex"    → has_g_quadruplex (boolean)
"uniprot_id"          → uniprot_id
"protein_name"        → protein_name
"gene_names"          → gene_names
"pubmed_id"           → pubmed_id
"buffer_conditions"   → buffer
"assay_description"   → assay_description
"similar_aptamers"    → similar_aptamers (array)
```

### Merge logic

Write a Python script `merge_databases.py` that:

1. Loads both datasets
2. Normalizes all records into a unified schema:

```json
{
  "id": "aptanexus_0001",
  "sequence": "AGCTCCAG...",
  "sequence_clean": "AGCTCCAG...",
  "aptamer_type": "DNA",
  "length": 56,
  "target_name": "β-conglutin",
  "target_type": "Protein",
  "kd_nM": 1.81,
  "pkd": 8.74,
  "gc_content": 0.536,
  "predicted_mfe": -26.22,
  "mfe_structure": ".(((((.(((......)).).)))))..",
  "has_g_quadruplex": false,
  "quality_tier": "P",
  "year": 2015,
  "doi": "10.1007/s00216-015-9179-z",
  "buffer": "10 mM phosphate, 138 mM NaCl...",
  "source": "aptanexus",
  "is_best": true
}
```

3. Computes derived features for ALL records:
   - `gc_content`: (G+C) / length
   - `a_freq`, `t_freq`, `g_freq`, `c_freq`: nucleotide frequencies
   - `purine_ratio`: (A+G) / length
   - `complexity`: Shannon entropy of dinucleotide frequencies
   - `has_g_quadruplex`: regex `G{3,}.{1,7}G{3,}.{1,7}G{3,}.{1,7}G{3,}`
   - `longest_repeat`: longest homopolymer run
   - `num_stems`: count of `(` in dot-bracket
   - `num_loops`: count of loop patterns `(\.+)` in dot-bracket

4. Deduplicates on (sequence, target_name):
   - If same sequence+target exists in both sources, keep the record with more fields filled
   - Mark `source: "both"` for duplicates
   - Expected overlap: ~500-800 records

5. Outputs `aptascope_merged.json` — the single data file the app imports

**Expected final count: ~13,000-13,500 unique records**

### How to run

```bash
python merge_databases.py \
  --aptanexus aptanexus_raw.jsonl \
  --aptadb aptadb_raw.json \
  --output src/data/aptascope_merged.json
```

Update the app to import from `aptascope_merged.json` instead of the current data file.

---

## BRANCH 1: Interaction Network Graph (Free)

### What it does

An interactive force-directed graph showing aptamer-target relationships. Nodes are targets (colored by type) and aptamers (colored by DNA/RNA). Edges represent binding interactions, with edge thickness proportional to affinity strength.

### Why it matters

Reveals patterns invisible in tables: which targets have the most aptamers, which aptamers cross-react with multiple targets, which target families cluster together, and where the field has gaps (orphan targets with only one known aptamer).

### Implementation

**Library:** D3.js v7 (import via CDN or npm)

**Data preparation:**
From the merged dataset, build a graph structure:

```javascript
{
  nodes: [
    { id: "target_thrombin", label: "Thrombin", type: "target", target_type: "Protein", aptamer_count: 47 },
    { id: "apt_GGTTGGTGTGGTTGG", label: "HD1", type: "aptamer", aptamer_type: "DNA", best_kd: 100 },
    ...
  ],
  links: [
    { source: "apt_GGTTGGTGTGGTTGG", target: "target_thrombin", kd_nM: 100, year: 1992 },
    ...
  ]
}
```

**Performance concern:** 13,000 records = 13,000 edges. D3 force simulation will choke. Solutions:
- Default view: show only targets with 5+ aptamers (~100 targets, ~2,000 edges)
- Filter controls: target type, affinity range, year range
- "Expand node" on click: show all aptamers for that target
- Use `d3-force` with `forceCollide` and `forceManyBody` with reduced charge for large graphs
- Use canvas rendering instead of SVG when node count > 500

**Interactions:**
- Hover node → show tooltip with name, type, count, best Kd
- Click target node → highlight all connected aptamers, show side panel with details
- Click aptamer node → show sequence, features, all targets it binds
- Zoom/pan with mouse
- Search box to find and center on a specific target
- Toggle: color by target_type / by affinity / by year
- Filter sidebar: target type checkboxes, aptamer count slider, affinity range

**Visual design:**
- Target nodes: circles, size proportional to log(aptamer_count), colored by target_type
- Aptamer nodes: small diamonds, blue for DNA, orange for RNA
- Edges: gray lines, opacity proportional to affinity strength (stronger = more opaque)
- Layout: force-directed with targets pulling their aptamers inward

**Add as new tab:** Explorer | Analytics | Analyzer | Targets | **Network**

---

## BRANCH 2: Secondary Structure Visualizer (Free)

### What it does

Renders aptamer secondary structures as interactive 2D diagrams from dot-bracket notation. 12,531 records already have this data — no computation needed.

### Why it matters

Researchers think about aptamers as shapes, not strings. A stem-loop diagram immediately communicates the binding architecture in a way dot-bracket notation doesn't. No existing aptamer database renders structures interactively in-browser.

### Implementation

**Option A — Use forna (recommended):**
The ViennaRNA team's `fornac` library renders RNA/DNA secondary structures from dot-bracket notation as interactive force-directed diagrams.

```bash
npm install fornac
```

```javascript
import { FornaContainer } from 'fornac';
// or load via CDN

const container = new FornaContainer('#structure-view', {
  animation: true,
  zoomable: true,
});
container.addRNA(dotBracketString, {
  sequence: aptamerSequence,
});
```

**Option B — Custom SVG rendering:**
If fornac is too heavy or hard to integrate with React, build a custom renderer:
1. Parse dot-bracket into paired/unpaired positions
2. Use Nussinov-style layout algorithm to place nucleotides in 2D
3. Render as SVG with:
   - Circles for nucleotides (colored A=green, T=red, G=amber, C=blue)
   - Lines for backbone
   - Arcs or lines for base pairs
   - Labels for position numbers

**Features:**
- Input: sequence + dot-bracket (auto-filled from database, or user-pasted)
- Render the 2D structure with nucleotide coloring
- Highlight regions on hover (stems, loops, bulges)
- Side-by-side comparison mode: user's sequence vs. best known binder
- Click a nucleotide → show position, whether it's paired, and with which position
- Export as SVG/PNG

**Integration points:**
- In the Analyzer view: after computing features, show the predicted structure
- In the Explorer: click any record → see its structure
- In Target Lookup: show structures of top binders side-by-side

---

## BRANCH 3: Aptamer Optimization Suggester (Free)

### What it does

Takes a user's sequence, identifies structural and compositional weaknesses, and suggests specific mutations based on statistical patterns from the 13,000-record dataset.

### Why it matters

This is the "so what do I do about it" feature. The Analyzer tells you your GC content is low — the Optimizer tells you which positions to change and why.

### Implementation

**Build a statistical model from the dataset (at build time, not runtime):**

1. For each target type (Protein, Small Molecule, Cell), compute:
   - Optimal GC content range (IQR of top-25% binders by Kd)
   - Optimal length range
   - Position-dependent nucleotide preferences (do high-affinity binders prefer G at position 1? C at the 3' end?)
   - Dinucleotide enrichment/depletion vs. random
   - Structural feature correlations (more stems = better? optimal loop size?)

2. Store as a precomputed JSON: `src/data/optimization_profiles.json`

**Runtime analysis of user's sequence:**

1. Identify weaknesses (deviations from optimal profiles):
   - GC content outside optimal range → flag
   - Homopolymer runs > 4nt → flag with position
   - Low complexity regions → flag
   - No predicted structure (all unpaired) → flag
   - Very negative or very positive MFE → flag

2. Generate suggestions with position-level specificity:
   ```
   Position 14-17 (AAAA): Homopolymer run reduces structural definition.
   → Suggest A16→G: introduces potential G-C pair without disrupting flanking stems.
   
   GC content 34% (below optimal 40-60% for protein targets):
   → Positions 8, 23, 31 are unpaired A/T in loop regions — 
      substituting any to G/C would raise GC without breaking predicted structure.
   
   3' tail (positions 38-45) is entirely unpaired:
   → Consider truncation to position 37 — 
      shorter aptamers with defined structure often bind tighter.
   ```

3. For each suggestion, check that it doesn't break predicted base pairs (use the dot-bracket notation to identify which positions are paired).

**UI:**
- Appears as a new panel in the Analyzer view, below the current results
- Each suggestion is a card with: position, current nt, suggested change, rationale, predicted impact
- "Apply suggestion" button → updates the sequence in the input and re-runs analysis
- "Apply all" → applies all non-conflicting suggestions at once

---

## BRANCH 4: Database Visualization Dashboard (Free)

### What it does

A rich, interactive overview of the entire 13,000-record dataset — the "state of aptamer science" view. Not per-aptamer analysis, but field-level patterns.

### Implementation

**Chart 1: Target landscape treemap**
Hierarchical treemap: top level = target type (Protein, Small Molecule, Cell...), second level = individual targets. Size = number of aptamers. Color = median affinity. Immediately shows where the field has invested attention.

**Chart 2: Affinity improvement over time**
Scatter plot: x = year, y = best Kd per target per year. Trend line showing whether aptamer affinities are improving over time (spoiler: they are, thanks to better SELEX methods). Overlay key method innovations (CE-SELEX, High-throughput SELEX) as annotations.

**Chart 3: Sequence space explorer**
Dimensionality reduction (t-SNE or UMAP) of all aptamer sequences based on k-mer frequency vectors. Compute at build time in Python, embed the 2D coordinates. Color by target type or affinity. Reveals whether protein-binders and small-molecule-binders occupy different regions of sequence space (they usually do).

For t-SNE computation (Python, at build time):
```python
from sklearn.manifold import TSNE
from sklearn.feature_extraction.text import CountVectorizer

# Extract 3-mer frequencies
def get_kmers(seq, k=3):
    return [seq[i:i+k] for i in range(len(seq)-k+1)]

kmer_strings = [' '.join(get_kmers(r['sequence'])) for r in records]
vectorizer = CountVectorizer()
X = vectorizer.fit_transform(kmer_strings)

tsne = TSNE(n_components=2, perplexity=30, random_state=42)
coords = tsne.fit_transform(X.toarray())

# Save coords alongside record IDs
```

**Chart 4: Quality tier breakdown**
Stacked bar chart by target type, showing A/B/C/P quality distribution. Helps researchers know how trustworthy the data is for their target class.

**Chart 5: Structure feature distributions**
Violin plots or histograms for MFE, stem count, loop count — split by high-affinity (<10 nM) vs. low-affinity (>1000 nM). Visually answers: "do good aptamers have more structure?"

**Add as new tab:** Explorer | Analytics | Analyzer | Targets | Network | **Dashboard**

---

## BRANCH 5: AI Sequence Analyzer (Requires API credits)

### What it does

User pastes a sequence → app computes features AND sends them to Claude → Claude writes a plain-English analysis interpreting the features in scientific context, comparing to known aptamers, and suggesting next experimental steps.

### Why it matters

Turns data into insight. A GC content of 52% means nothing to most researchers. "Your GC content is in the optimal range for protein binders, matching the profile of validated anti-VEGF aptamers" is actionable.

### Implementation

**Frontend:** Add a "Get AI Analysis" button in the Analyzer view, below the computed features.

**API call:**

```javascript
const response = await fetch("https://api.anthropic.com/v1/messages", {
  method: "POST",
  headers: { "Content-Type": "application/json" },
  body: JSON.stringify({
    model: "claude-sonnet-4-6",
    max_tokens: 1000,
    messages: [{
      role: "user",
      content: `You are an expert aptamer scientist. Analyze this candidate aptamer sequence and provide actionable insights.

SEQUENCE: ${sequence}
TYPE: ${dnaOrRna}
LENGTH: ${length} nt
GC CONTENT: ${gcContent}% (database median: ${medianGC}%)
PREDICTED MFE: ${mfe} kcal/mol
STRUCTURE: ${dotBracket}
G-QUADRUPLEX: ${hasGQuad}
COMPLEXITY (Shannon): ${entropy}

TOP 5 MOST SIMILAR KNOWN APTAMERS:
${similarAptamersTable}

TARGET (if specified by user): ${userTarget || "not specified"}

Provide:
1. Overall assessment (2-3 sentences)
2. Structural interpretation of the predicted fold
3. What the similarity matches suggest about likely target classes
4. Specific optimization suggestions with position numbers
5. Recommended next experimental steps

Be specific and quantitative. Reference the similar aptamers by name.`
    }]
  })
});
```

**Cost estimate:** ~$0.005 per analysis (Sonnet, ~500 input tokens + 500 output tokens). Budget-friendly for personal use.

**Streaming:** Use streaming response to show text appearing live — better UX than waiting.

**Fallback:** If no API key configured, show a banner: "AI analysis available — add your Anthropic API key in settings" and skip the call. The rest of the Analyzer works without it.

---

## BRANCH 6: AptaDiff Generative Design (Requires ML infrastructure)

### What it does

Integrates the AptaDiff diffusion model (https://github.com/wz-create/AptaDiff) to generate novel aptamer sequences conditioned on target properties.

### Why it matters

Moves the app from analysis to design — from "what exists" to "what should I make next." This is the frontier of computational aptamer science.

### Implementation path

**Step 1: Understand AptaDiff**
- Clone https://github.com/wz-create/AptaDiff
- Read the paper and code — it's a discrete diffusion model for DNA/RNA sequences
- It generates sequences conditioned on target embedding vectors
- Requires PyTorch

**Step 2: Build a Python API server**
- Wrap the AptaDiff model in a FastAPI server
- Endpoint: `POST /generate` — takes target name or embedding, returns N candidate sequences
- Run locally or deploy to a GPU cloud (Modal, Replicate, or RunPod)

```python
from fastapi import FastAPI
from aptadiff import AptaDiffModel

app = FastAPI()
model = AptaDiffModel.load_pretrained("checkpoints/best.pt")

@app.post("/generate")
def generate(target: str, n_candidates: int = 10):
    candidates = model.generate(target_embedding=encode_target(target), n=n_candidates)
    return {"candidates": [{"sequence": c.sequence, "confidence": c.score} for c in candidates]}
```

**Step 3: Connect to AptaScope frontend**
- New tab or section in Analyzer: "Generate Candidates"
- User selects a target → hits Generate → backend runs AptaDiff → returns candidates
- Each candidate is immediately analyzed with the existing Analyzer pipeline
- Display as a ranked table with all computed features + comparison to known binders

**Step 4: Validate loop**
- User picks promising candidates from the generated set
- App compares them against known binders: are the generated sequences novel? Do they have reasonable features?
- Export selected candidates for experimental testing

**Cost:** Free if running locally on a GPU laptop. ~$0.50-2.00 per generation batch on cloud GPU.

**Timeline:** This is a multi-day project. Don't attempt it tonight.

---

## Build order

```
Tonight (remaining hours):
  └── Phase 0: Database merge (30 min)
       └── Update app to use merged 13,000-record dataset

This week:
  ├── Branch 2: Structure Visualizer (free, ~3-4 hours)
  ├── Branch 3: Optimization Suggester (free, ~4-5 hours)
  └── Branch 1: Network Graph (free, ~5-6 hours)

Next week:
  └── Branch 4: Database Dashboard with t-SNE (free, ~4-5 hours)

When ready to invest:
  ├── Branch 5: AI Analyzer (~2-3 hours + API costs)
  └── Branch 6: AptaDiff integration (multi-day + GPU costs)
```

---

## For Claude Code

When starting a branch, tell Claude Code:

```
Read APTASCOPE_V2_SPEC.md, Phase 0 and Branch [N].
The merged data file is at src/data/aptascope_merged.json.
Build this branch as a new component and add it to the app navigation.
```

Each branch is independent — they can be built in any order and don't depend on each other (except Phase 0, which all branches need).
