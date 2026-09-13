// RNA/DNA secondary structure: dot-bracket parsing, a client-side Nussinov
// fold (for sequences the database didn't already fold with ViennaRNA —
// i.e. anything the user pastes into the Analyzer), and a 2D layout using
// the standard recursive circular-loop algorithm (the same technique tools
// like RNAplot/forna use): decompose the structure into a tree of loops,
// place each loop's nucleotides evenly around a circle sized to fit them,
// and extend each enclosed stem outward as a straight double-helix ladder
// toward its child loop. No pseudoknots are handled — dot-bracket notation
// (both ours and the database's ViennaRNA output) is pseudoknot-free by
// construction, so this is not a limitation in practice.

// --- dot-bracket <-> pair table -------------------------------------------

// pairTable[i] = index of the base i pairs with, or -1 if unpaired.
export function parseDotBracket(dotBracket) {
  const pairTable = new Array(dotBracket.length).fill(-1);
  const stack = [];
  for (let i = 0; i < dotBracket.length; i++) {
    const ch = dotBracket[i];
    if (ch === "(") {
      stack.push(i);
    } else if (ch === ")") {
      const j = stack.pop();
      if (j !== undefined) {
        pairTable[i] = j;
        pairTable[j] = i;
      }
    }
  }
  return pairTable;
}

// --- Nussinov fold (client-side fallback for arbitrary user sequences) ---

const MIN_HAIRPIN = 3; // minimum unpaired bases in a hairpin loop

function canPair(a, b) {
  const pairs = new Set(["AT", "TA", "AU", "UA", "GC", "CG", "GT", "TG", "GU", "UG"]);
  return pairs.has(a + b);
}

// Classic Nussinov maximum-base-pairing DP. Doesn't model real thermodynamics
// (no MFE energy value comes out of this — unlike the database's ViennaRNA-
// computed structures), but produces a valid, pseudoknot-free fold, which is
// exactly the fallback the original spec called for when a real folding
// engine isn't available.
export function nussinovFold(sequence) {
  const n = sequence.length;
  if (n < MIN_HAIRPIN + 2) return ".".repeat(n);

  const dp = Array.from({ length: n }, () => new Array(n).fill(0));

  for (let span = MIN_HAIRPIN + 1; span < n; span++) {
    for (let i = 0; i + span < n; i++) {
      const j = i + span;
      let best = dp[i + 1][j]; // i unpaired
      best = Math.max(best, dp[i][j - 1]); // j unpaired
      if (canPair(sequence[i], sequence[j])) {
        best = Math.max(best, dp[i + 1][j - 1] + 1);
      }
      for (let k = i + 1; k < j; k++) {
        best = Math.max(best, dp[i][k] + dp[k + 1][j]);
      }
      dp[i][j] = best;
    }
  }

  const pairTable = new Array(n).fill(-1);
  function traceback(i, j) {
    if (i >= j) return;
    if (dp[i][j] === dp[i + 1][j]) {
      traceback(i + 1, j);
    } else if (dp[i][j] === dp[i][j - 1]) {
      traceback(i, j - 1);
    } else if (canPair(sequence[i], sequence[j]) && dp[i][j] === dp[i + 1][j - 1] + 1) {
      pairTable[i] = j;
      pairTable[j] = i;
      traceback(i + 1, j - 1);
    } else {
      for (let k = i + 1; k < j; k++) {
        if (dp[i][j] === dp[i][k] + dp[k + 1][j]) {
          traceback(i, k);
          traceback(k + 1, j);
          return;
        }
      }
    }
  }
  traceback(0, n - 1);

  return pairTable.map((p, i) => (p === -1 ? "." : p > i ? "(" : ")")).join("");
}

// --- stems / loops (for hover highlighting) -------------------------------

// A stem is a maximal run of consecutive base pairs (i,j), (i+1,j-1), ...
// A loop is a maximal run of consecutive unpaired positions.
export function identifyRegions(pairTable) {
  const n = pairTable.length;
  const regionId = new Array(n).fill(-1);
  const visited = new Array(n).fill(false);
  const regions = [];

  let i = 0;
  while (i < n) {
    if (visited[i]) {
      i++;
      continue;
    }
    if (pairTable[i] === -1) {
      const start = i;
      while (i < n && pairTable[i] === -1) i++;
      const positions = range(start, i);
      const id = regions.length;
      regions.push({ type: "loop", positions });
      for (const p of positions) {
        regionId[p] = id;
        visited[p] = true;
      }
    } else if (pairTable[i] > i) {
      // Walk the 5' side of the stem forward while it keeps closing onto a
      // contiguous, still-open 3' partner (i.e. this is one continuous
      // double-helix ladder, not yet a bulge/loop/branch).
      const start = i;
      let j = pairTable[i];
      while (i + 1 < n && pairTable[i + 1] === j - 1 && i + 1 < pairTable[i + 1]) {
        i++;
        j = pairTable[i];
      }
      const positions = [];
      for (let k = start; k <= i; k++) positions.push(k);
      for (let k = pairTable[i]; k <= pairTable[start]; k++) positions.push(k);
      const id = regions.length;
      regions.push({ type: "stem", positions });
      for (const p of positions) {
        regionId[p] = id;
        visited[p] = true;
      }
      i++;
    } else {
      // Defensive fallback: a 3'-side partner reached before its 5' side
      // was marked visited shouldn't happen for a well-formed pair table,
      // but avoid corrupting state (or looping forever) if it ever does.
      regionId[i] = regions.length;
      regions.push({ type: "stem", positions: [i] });
      visited[i] = true;
      i++;
    }
  }
  return { regions, regionId };
}

function range(a, b) {
  const out = [];
  for (let k = a; k < b; k++) out.push(k);
  return out;
}

// --- recursive circular-loop layout ---------------------------------------

const POINT_SPACING = 15; // arc length (and stem rung spacing) per nucleotide — must exceed node diameter (2 * NODE_RADIUS in StructureViewer) or neighbors overlap
const RUNG_HALF_WIDTH = 7; // half the ladder width (5' strand <-> 3' strand)
const LOOP_GAP = 14; // gap between a stem's far end and the child loop's circle
const MIN_RADIUS = 16;

// Decompose a region of the backbone into an ordered list of "items" as seen
// walking around the loop that encloses it: each unpaired base is one item,
// and each stem (a maximal run of nested, contiguous base pairs) is one item
// whose enclosed region becomes a recursively-built child loop. This is the
// same loop/stem tree every RNA-structure-drawing algorithm is built on.
function buildLoopTree(pairTable, lo, hi) {
  const items = [];
  let i = lo;
  while (i <= hi) {
    if (pairTable[i] === -1) {
      items.push({ kind: "unpaired", pos: i });
      i++;
    } else {
      const outer5 = i;
      const outer3 = pairTable[i];
      let inner5 = i;
      while (
        inner5 + 1 <= hi &&
        pairTable[inner5 + 1] === pairTable[inner5] - 1 &&
        inner5 + 1 < pairTable[inner5 + 1]
      ) {
        inner5++;
      }
      const inner3 = pairTable[inner5];
      items.push({
        kind: "stem",
        outer5, outer3, inner5, inner3,
        child: buildLoopTree(pairTable, inner5 + 1, inner3 - 1),
      });
      i = outer3 + 1;
    }
  }
  return items;
}

// Circle radius needed to space `items` (plus one reserved slot for the
// entry point, if this loop has a parent stem) evenly at POINT_SPACING.
function loopRadius(items, hasParent) {
  const slots = items.length + (hasParent ? 1 : 0);
  return Math.max(MIN_RADIUS, (Math.max(slots, 3) * POINT_SPACING) / (2 * Math.PI));
}

function layoutLoop(items, positions, center, radius, baseAngle, hasParent) {
  const slots = items.length + (hasParent ? 1 : 0);
  const angleStep = (2 * Math.PI) / slots;
  // Slot 0 (at baseAngle) is the implicit entry/exit point where the parent
  // stem attaches — real items start one slot further round, so the loop
  // reads as a single unbroken arc away from the parent and back.
  const startIdx = hasParent ? 1 : 0;

  items.forEach((item, idx) => {
    const angle = baseAngle + (startIdx + idx) * angleStep;
    const cosA = Math.cos(angle);
    const sinA = Math.sin(angle);
    const circlePoint = { x: center.x + radius * cosA, y: center.y + radius * sinA };

    if (item.kind === "unpaired") {
      positions[item.pos] = circlePoint;
      return;
    }

    // Stem: place its outer pair straddling the circle point (along the
    // tangent), then extend both strands outward (along the radial
    // direction) as a straight ladder toward the child loop.
    const tangentX = -sinA, tangentY = cosA;
    const radialX = cosA, radialY = sinA;

    positions[item.outer5] = {
      x: circlePoint.x + tangentX * RUNG_HALF_WIDTH,
      y: circlePoint.y + tangentY * RUNG_HALF_WIDTH,
    };
    positions[item.outer3] = {
      x: circlePoint.x - tangentX * RUNG_HALF_WIDTH,
      y: circlePoint.y - tangentY * RUNG_HALF_WIDTH,
    };

    const numPairs = item.inner5 - item.outer5 + 1;
    let farCenter = circlePoint;
    for (let k = 1; k < numPairs; k++) {
      const t = k * POINT_SPACING;
      farCenter = { x: circlePoint.x + radialX * t, y: circlePoint.y + radialY * t };
      positions[item.outer5 + k] = {
        x: farCenter.x + tangentX * RUNG_HALF_WIDTH,
        y: farCenter.y + tangentY * RUNG_HALF_WIDTH,
      };
      positions[item.outer3 - k] = {
        x: farCenter.x - tangentX * RUNG_HALF_WIDTH,
        y: farCenter.y - tangentY * RUNG_HALF_WIDTH,
      };
    }

    if (item.child.length > 0) {
      const childRadius = loopRadius(item.child, true);
      const childCenter = {
        x: farCenter.x + radialX * (childRadius + LOOP_GAP),
        y: farCenter.y + radialY * (childRadius + LOOP_GAP),
      };
      // The child loop's entry slot (slot 0, at childBaseAngle) represents
      // the point on ITS OWN circle closest to this stem — i.e. it must
      // face back toward the parent (angle + PI), not away from it. Get
      // this backwards and the item diametrically opposite the entry (the
      // one pointing most "outward") ends up rotated by an extra 180° per
      // nesting level, which for chains of small loops/bulges can fold a
      // deeply-nested branch back on top of its own ancestor's geometry.
      layoutLoop(item.child, positions, childCenter, childRadius, angle + Math.PI, true);
    }
  });
}

export function computeLayout(sequence, pairTable) {
  const n = sequence.length;
  const positions = new Array(n);
  const tree = buildLoopTree(pairTable, 0, n - 1);
  const radius = loopRadius(tree, false);
  layoutLoop(tree, positions, { x: 0, y: 0 }, radius, -Math.PI / 2, false);
  return positions;
}
