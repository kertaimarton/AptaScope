// RCSB PDB lookup. None of the source databases carry a per-record PDB ID,
// so a target is matched to a solved structure by substring-checking its
// name against each PDB entry's title (e.g. target "VEGF" matches "G-rich
// VEGF aptamer with LNA modifications") — approximate, but the best link
// available without a real foreign key.
import pdbEnrichment from "../data/pdb_enrichment.json";

const MIN_MATCH_LEN = 4;

const ENTRIES = Object.values(pdbEnrichment).map((p) => ({
  ...p,
  _titleLower: p.title.toLowerCase(),
}));

export function findPdbStructures(targetName) {
  const lower = (targetName || "").toLowerCase();
  if (lower.length < MIN_MATCH_LEN) return [];
  return ENTRIES.filter((p) => p._titleLower.includes(lower));
}
