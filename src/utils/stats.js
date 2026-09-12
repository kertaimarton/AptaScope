// Percentile ranking and distribution helpers used by the Analyzer and
// Analytics views.

export function percentileRank(value, sortedValues) {
  if (sortedValues.length === 0) return null;
  let count = 0;
  for (const v of sortedValues) {
    if (v <= value) count += 1;
  }
  return (count / sortedValues.length) * 100;
}

export function median(values) {
  if (values.length === 0) return null;
  const sorted = [...values].sort((a, b) => a - b);
  const mid = Math.floor(sorted.length / 2);
  return sorted.length % 2 !== 0
    ? sorted[mid]
    : (sorted[mid - 1] + sorted[mid]) / 2;
}

export function quantile(sortedValues, q) {
  if (sortedValues.length === 0) return null;
  const pos = (sortedValues.length - 1) * q;
  const base = Math.floor(pos);
  const rest = pos - base;
  if (sortedValues[base + 1] !== undefined) {
    return sortedValues[base] + rest * (sortedValues[base + 1] - sortedValues[base]);
  }
  return sortedValues[base];
}

export function mean(values) {
  if (values.length === 0) return null;
  return values.reduce((sum, v) => sum + v, 0) / values.length;
}

export function ordinalSuffix(n) {
  const rounded = Math.round(n);
  const j = rounded % 10;
  const k = rounded % 100;
  if (j === 1 && k !== 11) return `${rounded}st`;
  if (j === 2 && k !== 12) return `${rounded}nd`;
  if (j === 3 && k !== 13) return `${rounded}rd`;
  return `${rounded}th`;
}

// Simple moving average over points already sorted by x, for trend lines.
export function movingAverage(points, windowSize = 15) {
  const sorted = [...points].sort((a, b) => a.x - b.x);
  const result = [];
  for (let i = 0; i < sorted.length; i++) {
    const start = Math.max(0, i - Math.floor(windowSize / 2));
    const end = Math.min(sorted.length, i + Math.ceil(windowSize / 2));
    const window = sorted.slice(start, end);
    const avgY = mean(window.map((p) => p.y));
    result.push({ x: sorted[i].x, y: avgY });
  }
  return result;
}
