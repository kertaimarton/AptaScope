// AlphaFold DB integration. AlphaFold predicts PROTEIN structures only — it
// has no model for nucleic acid folding, so this can never show an
// aptamer's own 3D structure. What it can show is the 3D structure of the
// protein TARGET an aptamer binds to, when that target carries a resolvable
// UniProt accession. The API is public, no key required, and CORS is open
// (access-control-allow-origin: *) on both the metadata endpoint and the
// structure file downloads, so this runs entirely client-side.

const API_BASE = "https://alphafold.ebi.ac.uk/api/prediction";

// A target's records may carry a UniProt ID two different ways depending on
// source (AptaDB's `uniprot_id`, or AptaNexus's `external_id` when
// `id_type` says it's a UniProt ID). Different records for the same target
// occasionally disagree (isoforms, inconsistent curation upstream) — take
// whichever accession appears most often among that target's records.
export function resolveUniProtId(records) {
  const counts = new Map();
  for (const r of records) {
    const uid = r.uniprot_id || (r.id_type === "UniProt ID" ? r.external_id : null);
    if (uid) counts.set(uid, (counts.get(uid) || 0) + 1);
  }
  if (counts.size === 0) return null;
  return [...counts.entries()].sort((a, b) => b[1] - a[1])[0][0];
}

export async function fetchAlphaFoldPrediction(uniprotId) {
  const res = await fetch(`${API_BASE}/${uniprotId}`);
  if (!res.ok) return null;
  const data = await res.json();
  return Array.isArray(data) && data.length > 0 ? data[0] : null;
}

// AlphaFold's own per-residue confidence (pLDDT) color scheme, as published
// alongside every AlphaFold DB structure — reusing it (rather than inventing
// a different one) means the coloring means the same thing here as it does
// on AlphaFold's own site.
export function plddtColor(plddt) {
  if (plddt > 90) return "#0053D6"; // very high
  if (plddt > 70) return "#65CBF3"; // confident
  if (plddt > 50) return "#FFDB13"; // low
  return "#FF7D45"; // very low
}
