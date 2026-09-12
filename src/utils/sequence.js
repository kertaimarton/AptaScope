// Sequence analysis functions — mirrors the feature computation in
// process_data.py so a user-pasted sequence is scored identically to the
// database records it's compared against.

const G_QUAD_PATTERN = /G{3,}.{1,7}G{3,}.{1,7}G{3,}.{1,7}G{3,}/;

export function detectType(sequence) {
  const hasU = sequence.includes("U");
  const hasT = sequence.includes("T");
  if (hasU && !hasT) return "RNA";
  return "DNA";
}

export function cleanSequence(raw) {
  return (raw || "")
    .toUpperCase()
    .replace(/[^ATGCU]/g, "");
}

export function isValidSequence(sequence) {
  return sequence.length > 0 && /^[ATGCU]+$/.test(sequence);
}

export function dinucleotideEntropy(seq) {
  if (seq.length < 2) return 0;
  const counts = new Map();
  for (let i = 0; i < seq.length - 1; i++) {
    const dinuc = seq.slice(i, i + 2);
    counts.set(dinuc, (counts.get(dinuc) || 0) + 1);
  }
  const total = seq.length - 1;
  let entropy = 0;
  for (const c of counts.values()) {
    const p = c / total;
    entropy -= p * Math.log2(p);
  }
  return entropy;
}

export function hasGQuadruplex(seq) {
  return G_QUAD_PATTERN.test(seq);
}

export function longestHomopolymerRun(seq) {
  if (!seq) return 0;
  let longest = 1;
  let current = 1;
  for (let i = 1; i < seq.length; i++) {
    if (seq[i] === seq[i - 1]) {
      current += 1;
      longest = Math.max(longest, current);
    } else {
      current = 1;
    }
  }
  return longest;
}

export function computeFeatures(sequence) {
  const type = detectType(sequence);
  const length = sequence.length;
  const counts = { A: 0, T: 0, U: 0, G: 0, C: 0 };
  for (const ch of sequence) {
    if (counts[ch] !== undefined) counts[ch] += 1;
  }
  const a = counts.A;
  const tOrU = type === "RNA" ? counts.U : counts.T;
  const g = counts.G;
  const c = counts.C;

  return {
    type,
    length,
    gc_content: (g + c) / length,
    a_freq: a / length,
    t_freq: tOrU / length,
    g_freq: g / length,
    c_freq: c / length,
    purine_ratio: (a + g) / length,
    complexity: dinucleotideEntropy(sequence),
    has_g_quadruplex: hasGQuadruplex(sequence),
    longest_repeat: longestHomopolymerRun(sequence),
  };
}

export function formatKd(kdNM) {
  if (kdNM === null || kdNM === undefined) return "—";
  if (kdNM < 1) return "< 1 nM";
  if (kdNM < 1000) return `${kdNM.toFixed(1)} nM`;
  if (kdNM < 1_000_000) return `${(kdNM / 1000).toFixed(2)} µM`;
  return `${(kdNM / 1_000_000).toFixed(2)} mM`;
}

export function kdColorClass(kdNM) {
  if (kdNM === null || kdNM === undefined) return "text-textsecondary";
  if (kdNM < 10) return "text-success";
  if (kdNM <= 1000) return "text-warning";
  return "text-danger";
}

const NUCLEOTIDE_COLORS = {
  A: "#238551",
  T: "#CD4246",
  U: "#CD4246",
  G: "#C87619",
  C: "#58A6FF",
};

export function nucleotideColor(base) {
  return NUCLEOTIDE_COLORS[base] || "#8F99A8";
}
