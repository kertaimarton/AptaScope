// Levenshtein edit distance and nearest-neighbor search against the
// database. O(n*m) per pair — pre-filter by length before scoring so we
// don't waste cycles comparing a 20nt input against a 150nt record.

export function levenshtein(a, b) {
  const m = a.length;
  const n = b.length;
  if (m === 0) return n;
  if (n === 0) return m;

  let prev = new Array(n + 1);
  let curr = new Array(n + 1);
  for (let j = 0; j <= n; j++) prev[j] = j;

  for (let i = 1; i <= m; i++) {
    curr[0] = i;
    const ai = a.charCodeAt(i - 1);
    for (let j = 1; j <= n; j++) {
      const cost = ai === b.charCodeAt(j - 1) ? 0 : 1;
      curr[j] = Math.min(
        prev[j] + 1, // deletion
        curr[j - 1] + 1, // insertion
        prev[j - 1] + cost // substitution
      );
    }
    [prev, curr] = [curr, prev];
  }
  return prev[n];
}

export function similarityPercent(distance, maxLen) {
  if (maxLen === 0) return 100;
  return Math.max(0, (1 - distance / maxLen) * 100);
}

// records: array of {sequence, target_name, kd_nM, source, id, ...}
export function findMostSimilar(inputSequence, records, topN = 5) {
  const inputLen = inputSequence.length;
  const minLen = inputLen * 0.5;
  const maxLen = inputLen * 2;

  const candidates = records.filter(
    (r) => r.sequence.length >= minLen && r.sequence.length <= maxLen
  );

  const scored = candidates.map((r) => {
    const distance = levenshtein(inputSequence, r.sequence);
    const longer = Math.max(inputLen, r.sequence.length);
    return {
      record: r,
      distance,
      similarity: similarityPercent(distance, longer),
    };
  });

  scored.sort((a, b) => b.similarity - a.similarity);
  return scored.slice(0, topN);
}
