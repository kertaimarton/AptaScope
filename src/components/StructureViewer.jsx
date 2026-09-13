import { useMemo, useRef, useState } from "react";
import { parseDotBracket, nussinovFold, identifyRegions, computeLayout } from "../utils/structure.js";
import { nucleotideColor } from "../utils/sequence.js";

const NODE_RADIUS = 6;

export default function StructureViewer({ sequence, dotBracket, label, height = 300 }) {
  const svgRef = useRef(null);
  const [hoverRegion, setHoverRegion] = useState(null);
  const [selected, setSelected] = useState(null);

  const { effectiveDotBracket, isPredicted } = useMemo(() => {
    if (dotBracket && dotBracket.length === sequence.length) {
      return { effectiveDotBracket: dotBracket, isPredicted: false };
    }
    return { effectiveDotBracket: nussinovFold(sequence), isPredicted: true };
  }, [sequence, dotBracket]);

  const pairTable = useMemo(() => parseDotBracket(effectiveDotBracket), [effectiveDotBracket]);
  const { regionId } = useMemo(() => identifyRegions(pairTable), [pairTable]);
  const layout = useMemo(
    () => computeLayout(sequence, pairTable),
    [sequence, pairTable]
  );

  const bounds = useMemo(() => {
    const xs = layout.map((p) => p.x);
    const ys = layout.map((p) => p.y);
    const pad = 16;
    const minX = Math.min(...xs) - pad;
    const minY = Math.min(...ys) - pad;
    const w = Math.max(...xs) - Math.min(...xs) + pad * 2;
    const h = Math.max(...ys) - Math.min(...ys) + pad * 2;
    return { minX, minY, w, h };
  }, [layout]);

  function exportSVG() {
    if (!svgRef.current) return;
    const serializer = new XMLSerializer();
    const source = serializer.serializeToString(svgRef.current);
    const blob = new Blob([source], { type: "image/svg+xml" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = "aptamer_structure.svg";
    a.click();
    URL.revokeObjectURL(url);
  }

  function exportPNG() {
    if (!svgRef.current) return;
    const serializer = new XMLSerializer();
    const source = serializer.serializeToString(svgRef.current);
    const svgBlob = new Blob([source], { type: "image/svg+xml" });
    const url = URL.createObjectURL(svgBlob);
    const img = new Image();
    img.onload = () => {
      const scale = 3;
      const canvas = document.createElement("canvas");
      canvas.width = bounds.w * scale;
      canvas.height = bounds.h * scale;
      const ctx = canvas.getContext("2d");
      ctx.fillStyle = "#111418";
      ctx.fillRect(0, 0, canvas.width, canvas.height);
      ctx.scale(scale, scale);
      ctx.drawImage(img, 0, 0, bounds.w, bounds.h);
      URL.revokeObjectURL(url);
      canvas.toBlob((blob) => {
        const pngUrl = URL.createObjectURL(blob);
        const a = document.createElement("a");
        a.href = pngUrl;
        a.download = "aptamer_structure.png";
        a.click();
        URL.revokeObjectURL(pngUrl);
      });
    };
    img.src = url;
  }

  const selectedInfo = selected !== null && selected !== undefined
    ? {
        position: selected + 1,
        base: sequence[selected],
        partner: pairTable[selected] !== -1 ? pairTable[selected] + 1 : null,
        partnerBase: pairTable[selected] !== -1 ? sequence[pairTable[selected]] : null,
      }
    : null;

  return (
    <div>
      {label && <div className="text-xs uppercase tracking-wider text-textsecondary mb-2">{label}</div>}
      <div style={{ width: "100%", height }}>
        <svg
          ref={svgRef}
          viewBox={`${bounds.minX} ${bounds.minY} ${bounds.w} ${bounds.h}`}
          width="100%"
          height="100%"
          preserveAspectRatio="xMidYMid meet"
          style={{ background: "#111418" }}
        >
          {/* backbone */}
          {layout.slice(0, -1).map((p, i) => {
            const q = layout[i + 1];
            const dim = hoverRegion !== null && regionId[i] !== hoverRegion && regionId[i + 1] !== hoverRegion;
            return (
              <line
                key={`bb-${i}`}
                x1={p.x} y1={p.y} x2={q.x} y2={q.y}
                stroke="#2F343C"
                strokeWidth={2}
                opacity={dim ? 0.15 : 1}
              />
            );
          })}
          {/* base pairs */}
          {pairTable.map((partner, i) => {
            if (partner <= i) return null;
            const p = layout[i];
            const q = layout[partner];
            const dim = hoverRegion !== null && regionId[i] !== hoverRegion;
            return (
              <line
                key={`bp-${i}`}
                x1={p.x} y1={p.y} x2={q.x} y2={q.y}
                stroke="#E0576B"
                strokeWidth={1.5}
                opacity={dim ? 0.1 : 0.85}
              />
            );
          })}
          {/* nucleotides */}
          {layout.map((p, i) => {
            const dim = hoverRegion !== null && regionId[i] !== hoverRegion;
            const isSelected = selected === i;
            return (
              <g
                key={`nt-${i}`}
                onMouseEnter={() => setHoverRegion(regionId[i])}
                onMouseLeave={() => setHoverRegion(null)}
                onClick={() => setSelected(i)}
                style={{ cursor: "pointer" }}
                opacity={dim ? 0.25 : 1}
              >
                <circle
                  cx={p.x} cy={p.y} r={NODE_RADIUS}
                  fill={nucleotideColor(sequence[i])}
                  stroke={isSelected ? "#F6F7F9" : "#111418"}
                  strokeWidth={isSelected ? 2 : 1}
                />
                <text
                  x={p.x} y={p.y}
                  textAnchor="middle"
                  dominantBaseline="central"
                  fontSize={7}
                  fontFamily="JetBrains Mono, monospace"
                  fill="#111418"
                  style={{ pointerEvents: "none" }}
                >
                  {sequence[i]}
                </text>
              </g>
            );
          })}
        </svg>
      </div>

      <div className="flex items-center justify-between mt-2 text-xs text-textsecondary">
        <span>
          {isPredicted
            ? "Predicted fold (base-pair maximization — no thermodynamic model in-browser)"
            : "Minimum free energy structure"}
          {selectedInfo && (
            <span className="ml-3 font-mono text-textprimary">
              pos {selectedInfo.position}: {selectedInfo.base}
              {selectedInfo.partner
                ? ` ↔ pos ${selectedInfo.partner} (${selectedInfo.partnerBase})`
                : " (unpaired)"}
            </span>
          )}
        </span>
        <span className="flex gap-2">
          <button
            onClick={exportSVG}
            className="px-2 py-0.5 border border-border hover:border-accent hover:text-textprimary"
          >
            Export SVG
          </button>
          <button
            onClick={exportPNG}
            className="px-2 py-0.5 border border-border hover:border-accent hover:text-textprimary"
          >
            Export PNG
          </button>
        </span>
      </div>
    </div>
  );
}
