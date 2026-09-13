import { findMostSimilar } from "./similarity.js";

self.onmessage = (e) => {
  const { sequence, records, topN } = e.data;
  const results = findMostSimilar(sequence, records, topN);
  // Strip the full record down to what the UI needs — smaller postMessage payload.
  const trimmed = results.map((r) => ({
    id: r.record.id,
    target_name: r.record.target_name,
    sequence: r.record.sequence,
    kd_nM: r.record.kd_nM,
    source: r.record.source,
    mfe_structure: r.record.mfe_structure,
    similarity: r.similarity,
    distance: r.distance,
  }));
  self.postMessage(trimmed);
};
