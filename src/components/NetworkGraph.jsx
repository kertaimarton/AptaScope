import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { forceSimulation, forceLink, forceManyBody, forceCollide, forceCenter, forceX, forceY } from "d3-force";
import { select } from "d3-selection";
import { zoom as d3zoom, zoomIdentity } from "d3-zoom";
import PageHeading from "./PageHeading.jsx";
import ColoredSequence, { CopyButton } from "./ColoredSequence.jsx";
import {
  buildTargetNodes,
  buildCrossLinks,
  findConnectedComponent,
  expandTarget,
  TARGET_TYPE_COLORS,
  kdToColor,
  yearToColor,
} from "../utils/networkGraph.js";
import { formatKd, kdColorClass } from "../utils/sequence.js";

const TARGET_TYPES = ["Protein", "Small Molecule", "Cell", "Nucleic Acid", "Microorganism", "Other"];
const COLOR_MODES = ["Target Type", "Affinity", "Year"];

function targetRadius(count) {
  return Math.max(5, Math.min(28, 4 + Math.log2(count + 1) * 3));
}

export default function NetworkGraph({ data }) {
  const containerRef = useRef(null);
  const canvasRef = useRef(null);
  const simRef = useRef(null);
  const nodesRef = useRef([]);
  const linksRef = useRef([]);
  const transformRef = useRef(zoomIdentity);
  const dimsRef = useRef({ width: 900, height: 640 });
  // The simulation's tick handler is registered once and must stay stable,
  // but `draw` is recreated whenever colorMode/hover/selected change (it
  // closes over them). Route through a ref so the tick handler always calls
  // the *current* draw — otherwise the simulation keeps invoking a stale
  // closure forever and things like the color-mode toggle silently never
  // take visual effect after the sim settles and stops ticking on its own.
  const drawRef = useRef(() => {});

  const [minCount, setMinCount] = useState(15);
  const [activeTypes, setActiveTypes] = useState(() => new Set(TARGET_TYPES));
  const [colorMode, setColorMode] = useState("Target Type");
  const [search, setSearch] = useState("");
  // null = overview (all filtered targets + cross-reactivity links between
  // them). Set = focus mode: only that one target + its full aptamer tree
  // is shown, everything else hidden — searching/clicking a target enters
  // focus mode; there's no "expand several at once" anymore, since that's
  // exactly what produced an unreadable scatter of disconnected clusters.
  const [focusedId, setFocusedId] = useState(null);
  const [selected, setSelected] = useState(null);
  const [hover, setHover] = useState(null); // { node, clientX, clientY }

  const allTargetNodes = useMemo(() => buildTargetNodes(data), [data]);
  const allTargetsById = useMemo(() => new Map(allTargetNodes.map((n) => [n.id, n])), [allTargetNodes]);
  const crossLinks = useMemo(() => buildCrossLinks(data), [data]);

  const yearBounds = useMemo(() => {
    const years = allTargetNodes.map((n) => n.year).filter(Boolean);
    return years.length ? { min: Math.min(...years), max: Math.max(...years) } : { min: 1990, max: 2025 };
  }, [allTargetNodes]);

  const filteredTargets = useMemo(
    () => allTargetNodes.filter((n) => n.count >= minCount && activeTypes.has(n.target_type)),
    [allTargetNodes, minCount, activeTypes]
  );

  const focusedTarget = focusedId ? allTargetsById.get(focusedId) : null;

  // Isolating just the clicked target would silently hide that it also
  // cross-reacts with others (e.g. Thrombin shares an aptamer with
  // Prothrombin, Gag-Pol polyprotein, and HIV Reverse Transcriptase in this
  // dataset) — so focus mode shows the whole connected group, not just the
  // one node clicked.
  const focusedCluster = useMemo(
    () => (focusedId ? findConnectedComponent(focusedId, crossLinks) : null),
    [focusedId, crossLinks]
  );

  const clusterTargets = useMemo(() => {
    if (!focusedCluster) return null;
    return [...focusedCluster].map((id) => allTargetsById.get(id)).filter(Boolean);
  }, [focusedCluster, allTargetsById]);

  // Memoized on the actual inputs (not just inline-computed) so this array
  // keeps a stable reference across renders that don't change focus or
  // filters — e.g. hover/select state updates — since it feeds a
  // useEffect below that restarts the force simulation whenever its
  // identity changes. Without this, hovering a node would re-trigger that
  // effect on every render and the layout would never settle.
  const displayTargets = useMemo(
    () => clusterTargets || filteredTargets,
    [clusterTargets, filteredTargets]
  );

  const nodeColor = useCallback(
    (n) => {
      if (n.kind === "aptamer") return n.aptamer_type === "DNA" ? "#58A6FF" : "#F0883E";
      if (colorMode === "Affinity") return kdToColor(n.medianKd);
      if (colorMode === "Year") return yearToColor(n.year, yearBounds.min, yearBounds.max);
      return TARGET_TYPE_COLORS[n.target_type] || "#8F99A8";
    },
    [colorMode, yearBounds]
  );

  // --- rebuild simulation nodes/links when filters or focus change ---
  useEffect(() => {
    const prevById = new Map(nodesRef.current.map((n) => [n.id, n]));
    const { width, height } = dimsRef.current;
    const newNodes = [];
    const newLinks = [];

    for (const t of displayTargets) {
      const prev = prevById.get(t.id);
      newNodes.push(
        prev
          ? Object.assign(prev, t)
          : { ...t, x: width / 2 + (Math.random() - 0.5) * 200, y: height / 2 + (Math.random() - 0.5) * 200 }
      );
    }

    if (clusterTargets) {
      // Expand every target in the connected component, not just the one
      // clicked — accumulating aptamerIdsSoFar across all of them means an
      // aptamer shared between e.g. Thrombin and Prothrombin becomes one
      // node with an edge to each, which is the actual cross-reactivity
      // this view exists to surface (not just "these are related" but
      // "here is the exact aptamer responsible").
      const aptamerIdsSoFar = new Set();
      for (const t of clusterTargets) {
        const { newAptamerNodes, newEdges } = expandTarget(t, aptamerIdsSoFar);
        const targetPos = prevById.get(t.id) || { x: width / 2, y: height / 2 };
        for (const an of newAptamerNodes) {
          aptamerIdsSoFar.add(an.id);
          const prev = prevById.get(an.id);
          newNodes.push(
            prev
              ? Object.assign(prev, an)
              : {
                  ...an,
                  x: targetPos.x + (Math.random() - 0.5) * 60,
                  y: targetPos.y + (Math.random() - 0.5) * 60,
                }
          );
        }
        newLinks.push(...newEdges.map((e) => ({ ...e, kind: "binding" })));
      }
    } else {
      const visibleIds = new Set(displayTargets.map((t) => t.id));
      for (const l of crossLinks) {
        if (visibleIds.has(l.source) && visibleIds.has(l.target)) {
          newLinks.push({ ...l, id: `${l.source}|${l.target}`, kind: "crosslink" });
        }
      }
    }

    nodesRef.current = newNodes;
    linksRef.current = newLinks;

    if (!simRef.current) {
      // forceCenter alone only recenters the average position — it does
      // nothing to stop mutual repulsion from spreading a few hundred nodes
      // far outside the canvas. forceX/forceY act as a weak spring back to
      // the middle that scales with distance, which is what actually keeps
      // the whole graph contained regardless of how many nodes are visible.
      simRef.current = forceSimulation(newNodes)
        .force("charge", forceManyBody().strength((d) => (d.kind === "target" ? -90 : -20)))
        .force("collide", forceCollide().radius((d) => (d.kind === "target" ? targetRadius(d.count) + 3 : 6)))
        .force("center", forceCenter(width / 2, height / 2))
        .force("x", forceX(width / 2).strength(0.05))
        .force("y", forceY(height / 2).strength(0.05))
        .on("tick", () => drawRef.current());
    } else {
      simRef.current.nodes(newNodes);
    }
    simRef.current.force(
      "link",
      forceLink(newLinks)
        .id((d) => d.id)
        .distance((l) => {
          if (l.kind === "crosslink") return 90;
          return l.target.kind === "target" ? targetRadius(l.target.count) + 30 : 60;
        })
        .strength((l) => (l.kind === "crosslink" ? 0.15 : 0.4))
    );
    simRef.current.alpha(0.6).restart();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [displayTargets, clusterTargets, crossLinks]);

  // --- canvas drawing ---
  const draw = useCallback(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const ctx = canvas.getContext("2d");
    const { width, height } = dimsRef.current;
    const dpr = window.devicePixelRatio || 1;
    ctx.save();
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.fillStyle = "#111418";
    ctx.fillRect(0, 0, width, height);

    const t = transformRef.current;
    ctx.translate(t.x, t.y);
    ctx.scale(t.k, t.k);

    // edges
    for (const l of linksRef.current) {
      const s = typeof l.source === "object" ? l.source : nodesRef.current.find((n) => n.id === l.source);
      const e = typeof l.target === "object" ? l.target : nodesRef.current.find((n) => n.id === l.target);
      if (!s || !e || s.x === undefined || e.x === undefined) continue;
      ctx.beginPath();
      ctx.moveTo(s.x, s.y);
      ctx.lineTo(e.x, e.y);
      if (l.kind === "crosslink") {
        ctx.strokeStyle = "rgba(224,87,107,0.5)";
        ctx.lineWidth = Math.min(3, 1 + l.sharedCount * 0.4) / t.k;
        ctx.setLineDash([4 / t.k, 3 / t.k]);
      } else {
        const opacity = l.kd_nM ? Math.max(0.15, Math.min(0.7, 1 - Math.log10(l.kd_nM + 1) / 6)) : 0.15;
        ctx.strokeStyle = `rgba(143,153,168,${opacity})`;
        ctx.lineWidth = 1 / t.k;
        ctx.setLineDash([]);
      }
      ctx.stroke();
    }
    ctx.setLineDash([]);

    // nodes
    for (const n of nodesRef.current) {
      if (n.x === undefined) continue;
      const color = nodeColor(n);
      const isHover = hover?.node?.id === n.id;
      const isSelected = selected?.id === n.id;
      ctx.fillStyle = color;
      ctx.strokeStyle = isSelected ? "#F6F7F9" : isHover ? "#F6F7F9" : "#111418";
      ctx.lineWidth = (isSelected ? 2.5 : isHover ? 2 : 1) / t.k;

      if (n.kind === "target") {
        const r = targetRadius(n.count);
        ctx.beginPath();
        ctx.arc(n.x, n.y, r, 0, 2 * Math.PI);
        ctx.fill();
        ctx.stroke();
      } else {
        const r = 5;
        ctx.beginPath();
        ctx.moveTo(n.x, n.y - r);
        ctx.lineTo(n.x + r, n.y);
        ctx.lineTo(n.x, n.y + r);
        ctx.lineTo(n.x - r, n.y);
        ctx.closePath();
        ctx.fill();
        ctx.stroke();
      }
    }
    ctx.restore();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [nodeColor, hover, selected]);

  useEffect(() => {
    drawRef.current = draw;
    draw();
  }, [draw]);

  // --- zoom/pan ---
  useEffect(() => {
    const canvas = canvasRef.current;
    const behavior = d3zoom()
      .scaleExtent([0.15, 6])
      .on("zoom", (event) => {
        transformRef.current = event.transform;
        drawRef.current();
      });
    select(canvas).call(behavior);
    canvas.__zoomBehavior = behavior;
    return () => select(canvas).on(".zoom", null);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  // --- responsive sizing ---
  // A one-shot resize on mount can fire before the surrounding CSS grid has
  // settled its column widths, measuring a stale (too-small) clientWidth
  // and leaving the canvas short of its container — ResizeObserver instead
  // reports the container's actual size whenever it changes, including the
  // first real layout pass.
  useEffect(() => {
    function resize(width) {
      if (!canvasRef.current) return;
      const height = 640;
      dimsRef.current = { width, height };
      const dpr = window.devicePixelRatio || 1;
      canvasRef.current.width = width * dpr;
      canvasRef.current.height = height * dpr;
      canvasRef.current.style.width = `${width}px`;
      canvasRef.current.style.height = `${height}px`;
      simRef.current?.force("center", forceCenter(width / 2, height / 2));
      simRef.current?.force("x", forceX(width / 2).strength(0.05));
      simRef.current?.force("y", forceY(height / 2).strength(0.05));
      drawRef.current();
    }
    const observer = new ResizeObserver((entries) => {
      const width = entries[0]?.contentRect.width;
      if (width) resize(width);
    });
    if (containerRef.current) observer.observe(containerRef.current);
    return () => observer.disconnect();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  // --- pointer interactions: hover + click ---
  function nodeAtPoint(clientX, clientY) {
    const canvas = canvasRef.current;
    const rect = canvas.getBoundingClientRect();
    const t = transformRef.current;
    const x = (clientX - rect.left - t.x) / t.k;
    const y = (clientY - rect.top - t.y) / t.k;
    let best = null;
    let bestDist = Infinity;
    for (const n of nodesRef.current) {
      if (n.x === undefined) continue;
      const r = n.kind === "target" ? targetRadius(n.count) : 6;
      const dist = Math.hypot(n.x - x, n.y - y);
      if (dist <= r + 2 && dist < bestDist) {
        best = n;
        bestDist = dist;
      }
    }
    return best;
  }

  function handleMouseMove(e) {
    const node = nodeAtPoint(e.clientX, e.clientY);
    setHover(node ? { node, clientX: e.clientX, clientY: e.clientY } : null);
  }

  function handleClick(e) {
    const node = nodeAtPoint(e.clientX, e.clientY);
    if (!node) {
      // Clicking empty canvas backs out of focus mode, back to the overview.
      setSelected(null);
      setFocusedId(null);
      return;
    }
    setSelected(node);
    if (node.kind === "target") {
      setFocusedId((prev) => (prev === node.id ? null : node.id));
    }
  }

  function backToOverview() {
    setFocusedId(null);
    setSelected(null);
  }

  function centerOn(node) {
    const { width, height } = dimsRef.current;
    const k = 1.5;
    const newTransform = zoomIdentity.translate(width / 2 - node.x * k, height / 2 - node.y * k).scale(k);
    select(canvasRef.current)
      .transition()
      .duration(500)
      .call(canvasRef.current.__zoomBehavior.transform, newTransform);
  }

  function handleSearchSelect(targetNode) {
    setSelected(targetNode);
    setFocusedId(targetNode.id);
    // wait a tick for the simulation to place the node before centering
    setTimeout(() => centerOn(nodesRef.current.find((n) => n.id === targetNode.id) || targetNode), 300);
  }

  const searchResults = useMemo(() => {
    const q = search.trim().toLowerCase();
    if (!q) return [];
    return allTargetNodes.filter((n) => n.label.toLowerCase().includes(q)).slice(0, 8);
  }, [search, allTargetNodes]);

  function toggleType(type) {
    setActiveTypes((prev) => {
      const next = new Set(prev);
      if (next.has(type)) next.delete(type);
      else next.add(type);
      return next;
    });
  }

  return (
    <div className="space-y-4">
      <PageHeading>Network</PageHeading>
      <p className="text-[13px] text-textsecondary">
        {clusterTargets ? (
          <>
            Showing every aptamer bound to{" "}
            <span className="text-textprimary">{focusedTarget.label}</span>
            {clusterTargets.length > 1 ? (
              <>
                {" "}
                and its cross-reactive group —{" "}
                <span className="text-textprimary">
                  {clusterTargets
                    .filter((t) => t.id !== focusedTarget.id)
                    .map((t) => t.label)
                    .join(", ")}
                </span>{" "}
                share at least one aptamer with it.
              </>
            ) : (
              "."
            )}{" "}
            Click empty space, the center node, or "Back to overview" to return.
          </>
        ) : (
          <>
            Every target is shown, sized by aptamer count; dashed rose lines connect targets
            that share a cross-reactive aptamer. Click (or search for) a target to focus on
            its full aptamer tree. Scroll to zoom, drag to pan.
          </>
        )}
      </p>

      <div className="grid grid-cols-1 lg:grid-cols-4 gap-4">
        {/* Filter sidebar */}
        <div className="bg-surface border border-border rounded-md p-4 space-y-4 lg:col-span-1 h-fit">
          <div>
            <label className="block text-xs uppercase tracking-wider text-textsecondary mb-1">
              Search targets
            </label>
            <input
              type="text"
              value={search}
              onChange={(e) => setSearch(e.target.value)}
              placeholder="e.g. thrombin"
              className="w-full bg-bg border border-border rounded px-3 py-1.5 text-sm focus:outline-none focus:border-accent"
            />
            {searchResults.length > 0 && (
              <div className="mt-1 border border-border rounded overflow-hidden">
                {searchResults.map((n) => (
                  <button
                    key={n.id}
                    onClick={() => {
                      setSearch("");
                      handleSearchSelect(n);
                    }}
                    className="block w-full text-left px-2 py-1 text-xs hover:bg-bg"
                  >
                    {n.label} <span className="text-textsecondary font-mono">({n.count})</span>
                  </button>
                ))}
              </div>
            )}
          </div>

          {focusedTarget ? (
            <button
              onClick={backToOverview}
              className="text-xs uppercase tracking-wider border border-accent text-accent px-3 py-1.5 w-full hover:bg-accent hover:text-bg transition-colors"
            >
              ← Back to overview
            </button>
          ) : (
            <>
              <div>
                <label className="block text-xs uppercase tracking-wider text-textsecondary mb-1">
                  Min aptamers per target: <span className="font-mono text-textprimary">{minCount}</span>
                </label>
                <input
                  type="range"
                  min={1}
                  max={100}
                  value={minCount}
                  onChange={(e) => setMinCount(Number(e.target.value))}
                  className="w-full accent-accent"
                />
                <div className="text-xs text-textsecondary font-mono">
                  {filteredTargets.length} targets shown
                </div>
              </div>

              <div>
                <label className="block text-xs uppercase tracking-wider text-textsecondary mb-2">
                  Target type
                </label>
                <div className="space-y-1">
                  {TARGET_TYPES.map((type) => (
                    <label key={type} className="flex items-center gap-2 text-sm cursor-pointer">
                      <input
                        type="checkbox"
                        checked={activeTypes.has(type)}
                        onChange={() => toggleType(type)}
                        className="accent-accent"
                      />
                      <span
                        className="inline-block w-2.5 h-2.5 shrink-0"
                        style={{ backgroundColor: TARGET_TYPE_COLORS[type] }}
                      />
                      {type}
                    </label>
                  ))}
                </div>
              </div>
            </>
          )}

          <div>
            <label className="block text-xs uppercase tracking-wider text-textsecondary mb-2">
              Color targets by
            </label>
            <div className="flex flex-col gap-1">
              {COLOR_MODES.map((mode) => (
                <button
                  key={mode}
                  onClick={() => setColorMode(mode)}
                  className={`text-left text-sm px-2 py-1 rounded border ${
                    colorMode === mode
                      ? "border-accent text-accent"
                      : "border-border text-textsecondary hover:text-textprimary"
                  }`}
                >
                  {mode}
                </button>
              ))}
            </div>
          </div>
        </div>

        {/* Graph canvas */}
        <div className="lg:col-span-3 space-y-4">
          <div
            ref={containerRef}
            className="bg-surface border border-border rounded-md relative overflow-hidden"
          >
            <canvas
              ref={canvasRef}
              onMouseMove={handleMouseMove}
              onMouseLeave={() => setHover(null)}
              onClick={handleClick}
              style={{ display: "block", cursor: "pointer" }}
            />
            {hover && (
              <div
                className="fixed pointer-events-none bg-bg border border-border rounded px-2 py-1.5 text-xs z-50 shadow-lg"
                style={{ left: hover.clientX + 12, top: hover.clientY + 12 }}
              >
                {hover.node.kind === "target" ? (
                  <>
                    <div className="font-medium text-textprimary">{hover.node.label}</div>
                    <div className="text-textsecondary">{hover.node.target_type}</div>
                    <div className="text-textsecondary">
                      {hover.node.count} aptamers · best{" "}
                      <span className={kdColorClass(hover.node.bestKd)}>{formatKd(hover.node.bestKd)}</span>
                    </div>
                  </>
                ) : (
                  <>
                    <div className="font-mono text-textprimary">
                      {hover.node.sequence.slice(0, 24)}
                      {hover.node.sequence.length > 24 ? "…" : ""}
                    </div>
                    <div className={kdColorClass(hover.node.bestKd)}>{formatKd(hover.node.bestKd)}</div>
                  </>
                )}
              </div>
            )}
          </div>

          {/* Side panel */}
          {selected && (
            <div className="bg-surface border border-border rounded-md p-4">
              {selected.kind === "target" ? <TargetPanel node={selected} /> : <AptamerPanel node={selected} />}
            </div>
          )}
        </div>
      </div>
    </div>
  );
}

function TargetPanel({ node }) {
  return (
    <div>
      <div className="flex items-center justify-between mb-2">
        <h3 className="text-sm font-semibold">{node.label}</h3>
        <span
          className="text-xs px-2 py-0.5 rounded-full border"
          style={{ borderColor: TARGET_TYPE_COLORS[node.target_type], color: TARGET_TYPE_COLORS[node.target_type] }}
        >
          {node.target_type}
        </span>
      </div>
      <div className="grid grid-cols-3 gap-4 text-sm">
        <div>
          <div className="text-xs uppercase tracking-wider text-textsecondary">Aptamers</div>
          <div className="font-mono">{node.count}</div>
        </div>
        <div>
          <div className="text-xs uppercase tracking-wider text-textsecondary">Best Kd</div>
          <div className={`font-mono ${kdColorClass(node.bestKd)}`}>{formatKd(node.bestKd)}</div>
        </div>
        <div>
          <div className="text-xs uppercase tracking-wider text-textsecondary">Median Year</div>
          <div className="font-mono">{node.year || "—"}</div>
        </div>
      </div>
    </div>
  );
}

function AptamerPanel({ node }) {
  const r = node.record;
  return (
    <div className="space-y-2">
      <div className="flex items-center justify-between">
        <h3 className="text-sm font-semibold">Aptamer</h3>
        <span
          className={`text-xs px-2 py-0.5 rounded-full border ${
            node.aptamer_type === "DNA" ? "border-dna text-dna" : "border-rna text-rna"
          }`}
        >
          {node.aptamer_type}
        </span>
      </div>
      <div className="flex items-center gap-2">
        <ColoredSequence sequence={node.sequence} className="text-sm" />
        <CopyButton text={node.sequence} />
      </div>
      <div className="grid grid-cols-3 gap-4 text-sm">
        <div>
          <div className="text-xs uppercase tracking-wider text-textsecondary">Length</div>
          <div className="font-mono">{r.length} nt</div>
        </div>
        <div>
          <div className="text-xs uppercase tracking-wider text-textsecondary">Kd</div>
          <div className={`font-mono ${kdColorClass(r.kd_nM)}`}>{formatKd(r.kd_nM)}</div>
        </div>
        <div>
          <div className="text-xs uppercase tracking-wider text-textsecondary">Target</div>
          <div>{r.target_name}</div>
        </div>
      </div>
    </div>
  );
}
