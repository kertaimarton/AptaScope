// Branch 3 — Aptamer Optimization Suggester.
//
// Compares a user's sequence against the statistical profile of top-25%
// binders for a given target type (src/data/optimization_profiles.json,
// built by build_optimization_profiles.py from the 13,306-record dataset)
// and produces position-specific mutation suggestions, each carrying a
// precomputed `newSequence` so the UI can apply it with a single click.
// Every suggested substitution is checked against the predicted pair table
// first — a position that's base-paired only gets a suggestion if the
// replacement can still pair with its existing partner, per the spec's
// "check that it doesn't break predicted base pairs" requirement.

import { canPair } from "./structure.js";

const NUM_BINS = 5;
const BASES = ["A", "T", "G", "C"];
const MIN_HOMOPOLYMER = 5;
const MIN_TAIL = 6;

function positionBin(i, length) {
  const frac = i / length;
  return Math.min(NUM_BINS - 1, Math.floor(frac * NUM_BINS));
}

function norm(ch) {
  return ch === "U" ? "T" : ch;
}

// Best replacement base for `pos` per the profile's position-preference
// enrichment, excluding the current base, that still preserves pairing if
// `pos` is base-paired. Returns null if no valid alternative exists.
function bestReplacement(sequence, pos, pairTable, profile) {
  const bin = positionBin(pos, sequence.length);
  const enrichment = profile?.position_preferences?.[bin]?.enrichment;
  if (!enrichment) return null;
  const current = norm(sequence[pos]);
  const isRNA = sequence.includes("U");
  const candidates = BASES.filter((b) => b !== current)
    .map((b) => ({ base: isRNA && b === "T" ? "U" : b, score: enrichment[b] }))
    .filter((c) => c.score !== null && c.score !== undefined)
    .sort((a, b) => b.score - a.score);

  const partner = pairTable[pos];
  for (const c of candidates) {
    if (partner === -1 || partner === undefined) return c;
    if (canPair(c.base, sequence[partner])) return c;
  }
  return null;
}

function findHomopolymerRuns(sequence) {
  const runs = [];
  let start = 0;
  for (let i = 1; i <= sequence.length; i++) {
    if (i === sequence.length || sequence[i] !== sequence[start]) {
      if (i - start >= MIN_HOMOPOLYMER) runs.push({ start, end: i - 1, base: sequence[start] });
      start = i;
    }
  }
  return runs.sort((a, b) => b.end - b.start - (a.end - a.start));
}

function homopolymerSuggestions(sequence, pairTable, profile) {
  const out = [];
  for (const run of findHomopolymerRuns(sequence).slice(0, 3)) {
    const pos = run.start + Math.floor((run.end - run.start) / 2);
    const replacement = bestReplacement(sequence, pos, pairTable, profile);
    const runLabel = `${pos1(run.start)}-${pos1(run.end)}`;
    if (!replacement) {
      out.push({
        id: `homopolymer-${run.start}`,
        category: "homopolymer",
        title: `Homopolymer run at ${runLabel} (${run.base.repeat(run.end - run.start + 1)})`,
        rationale:
          "Runs longer than 4nt reduce structural definition and can promote misfolding, but every position in this run is base-paired with no substitution that preserves the pair — a manual redesign of this stem is needed.",
        impact: "Flag only",
        positions: [pos],
        currentFragment: null,
        suggestedFragment: null,
        newSequence: null,
      });
      continue;
    }
    const newSequence = sequence.slice(0, pos) + replacement.base + sequence.slice(pos + 1);
    out.push({
      id: `homopolymer-${run.start}`,
      category: "homopolymer",
      title: `Homopolymer run at ${runLabel} (${run.base.repeat(run.end - run.start + 1)})`,
      rationale: `Runs longer than 4nt reduce structural definition. Position ${pos1(pos)} (${sequence[pos]}) is ${
        pairTable[pos] === -1 ? "unpaired" : "base-paired"
      } and top binders favor ${replacement.base} there (${replacement.score.toFixed(2)}x enrichment vs. the rest of the cohort).`,
      impact: `Breaks up the run without disrupting ${pairTable[pos] === -1 ? "structure" : "the existing base pair"}`,
      positions: [pos],
      currentFragment: sequence[pos],
      suggestedFragment: replacement.base,
      newSequence,
    });
  }
  return out;
}

function gcSuggestion(sequence, features, pairTable, profile) {
  if (!profile?.gc_range) return null;
  const [lo, hi] = profile.gc_range;
  if (features.gc_content >= lo && features.gc_content <= hi) return null;

  const tooLow = features.gc_content < lo;
  const targetCount = Math.ceil((tooLow ? lo : hi) * sequence.length);
  const currentCount = Math.round(features.gc_content * sequence.length);
  const needed = Math.abs(targetCount - currentCount);
  if (needed === 0) return null;

  // Candidate unpaired positions currently holding the "wrong" kind of base
  // (A/T when we need more GC, G/C when we need less), ranked by how
  // strongly top binders favor a swap at that position's bin.
  const candidates = [];
  for (let i = 0; i < sequence.length; i++) {
    if (pairTable[i] !== -1) continue; // don't touch paired positions
    const isGC = sequence[i] === "G" || sequence[i] === "C";
    if (tooLow && isGC) continue;
    if (!tooLow && !isGC) continue;
    const replacement = bestReplacement(sequence, i, pairTable, profile);
    if (!replacement) continue;
    const replacementIsGC = replacement.base === "G" || replacement.base === "C";
    if (tooLow && !replacementIsGC) continue;
    if (!tooLow && replacementIsGC) continue;
    candidates.push({ pos: i, base: replacement.base, score: replacement.score });
  }
  candidates.sort((a, b) => b.score - a.score);
  const chosen = candidates.slice(0, Math.min(needed, 4));
  if (chosen.length === 0) return null;

  let seqArr = sequence.split("");
  for (const c of chosen) seqArr[c.pos] = c.base;
  const newSequence = seqArr.join("");
  const positions = chosen.map((c) => c.pos);
  const newGC = (newSequence.match(/[GC]/g) || []).length / newSequence.length;

  return {
    id: "gc-content",
    category: "gc",
    title: `GC content ${(features.gc_content * 100).toFixed(0)}% is ${tooLow ? "below" : "above"} the ${(lo * 100).toFixed(
      0
    )}–${(hi * 100).toFixed(0)}% range typical of this target type's top binders`,
    rationale: `Positions ${positions.map(pos1).join(", ")} are unpaired and ${
      tooLow ? "A/T" : "G/C"
    } in loop/bulge regions — substituting ${tooLow ? "them to G/C" : "them to A/T"} moves GC content toward the optimal range without breaking predicted structure.`,
    impact: `GC ${(features.gc_content * 100).toFixed(0)}% → ${(newGC * 100).toFixed(0)}%`,
    positions,
    currentFragment: positions.map((p) => sequence[p]).join(""),
    suggestedFragment: chosen.map((c) => c.base).join(""),
    newSequence,
  };
}

function tailSuggestion(sequence, pairTable) {
  let tailStart = sequence.length;
  for (let i = sequence.length - 1; i >= 0; i--) {
    if (pairTable[i] !== -1) break;
    tailStart = i;
  }
  const tailLen = sequence.length - tailStart;
  if (tailLen < MIN_TAIL) return null;

  const newSequence = sequence.slice(0, tailStart);
  return {
    id: "tail-truncation",
    category: "tail",
    title: `3' tail (positions ${pos1(tailStart)}-${pos1(sequence.length - 1)}) is entirely unpaired`,
    rationale:
      "Shorter aptamers with fully defined structure often bind tighter than the same sequence with a floppy, unstructured tail — consider truncating to the last paired position.",
    impact: `${sequence.length}nt → ${newSequence.length}nt`,
    positions: null,
    currentFragment: sequence.slice(tailStart),
    suggestedFragment: "(removed)",
    newSequence,
  };
}

function structureFlag(pairTable) {
  const pairedCount = pairTable.filter((p) => p !== -1).length;
  if (pairedCount > 0) return null;
  return {
    id: "no-structure",
    category: "structure",
    title: "No predicted secondary structure",
    rationale:
      "This sequence folds as entirely unpaired under the fallback predictor. Most validated aptamers rely on a stem-loop or G-quadruplex fold to present a rigid binding surface — a fully flexible sequence is unusual among strong binders.",
    impact: "Flag only",
    positions: null,
    currentFragment: null,
    suggestedFragment: null,
    newSequence: null,
  };
}

function pos1(i) {
  return i + 1; // 1-indexed for display
}

export function generateOptimizationSuggestions(sequence, features, pairTable, profile) {
  const suggestions = [];
  const structFlag = structureFlag(pairTable);
  if (structFlag) suggestions.push(structFlag);

  const gc = gcSuggestion(sequence, features, pairTable, profile);
  if (gc) suggestions.push(gc);

  suggestions.push(...homopolymerSuggestions(sequence, pairTable, profile));

  const tail = tailSuggestion(sequence, pairTable);
  if (tail) suggestions.push(tail);

  return suggestions;
}

// Greedily combine every substitution-type suggestion whose positions don't
// overlap an already-applied one (first in list wins), then apply at most
// one truncation last — it only removes trailing bases, so it can't
// invalidate earlier positions' edits.
export function applyAllSuggestions(sequence, suggestions) {
  const used = new Set();
  let seqArr = sequence.split("");
  let truncateAt = null;

  for (const s of suggestions) {
    if (!s.newSequence) continue;
    if (s.category === "tail") {
      truncateAt = s.newSequence.length;
      continue;
    }
    if (!s.positions) continue;
    if (s.positions.some((p) => used.has(p))) continue;
    const suggestedBases = s.suggestedFragment;
    s.positions.forEach((p, idx) => {
      seqArr[p] = suggestedBases[idx];
      used.add(p);
    });
  }

  let result = seqArr.join("");
  if (truncateAt !== null) result = result.slice(0, truncateAt);
  return result;
}
