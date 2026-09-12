import { useEffect, useMemo, useRef, useState } from "react";
import {
  ResponsiveContainer,
  RadarChart,
  PolarGrid,
  PolarAngleAxis,
  PolarRadiusAxis,
  Radar,
  Tooltip,
  Legend,
} from "recharts";
import ColoredSequence, { CopyButton } from "./ColoredSequence.jsx";
import { cleanSequence, computeFeatures, formatKd, kdColorClass } from "../utils/sequence.js";
import { percentileRank, median, ordinalSuffix } from "../utils/stats.js";
import PageHeading from "./PageHeading.jsx";

const EXAMPLES = [
  { label: "Thrombin Aptamer (HD1)", sequence: "GGTTGGTGTGGTTGG" },
  { label: "VEGF Aptamer", sequence: "TGTGGGGGTGGACGGGCCGGGTAGA" },
];

function buildDistributions(data) {
  const dist = {
    gc_content: [],
    complexity: [],
    length: [],
    purine_ratio: [],
    predicted_mfe: [],
  };
  for (const r of data) {
    dist.gc_content.push(r.gc_content);
    dist.complexity.push(r.complexity);
    dist.length.push(r.length);
    dist.purine_ratio.push(r.purine_ratio);
    if (r.predicted_mfe !== null) dist.predicted_mfe.push(r.predicted_mfe);
  }
  for (const key of Object.keys(dist)) dist[key].sort((a, b) => a - b);
  return dist;
}

function PercentileGauge({ percentile }) {
  const pct = percentile === null ? 0 : Math.max(0, Math.min(100, percentile));
  return (
    <div className="h-1.5 bg-bg rounded-full overflow-hidden border border-border">
      <div
        className="h-full bg-accent"
        style={{ width: `${pct}%` }}
      />
    </div>
  );
}

function generateSuggestions(features) {
  const suggestions = [];
  if (features.gc_content < 0.35 || features.gc_content > 0.7) {
    suggestions.push(
      "GC content is outside the typical range for high-affinity binders (40–60%). Consider adjusting."
    );
  }
  if (features.predicted_mfe === null || features.predicted_mfe >= -2) {
    suggestions.push(
      "Low structural stability. High-affinity aptamers typically form defined secondary structures."
    );
  }
  if (features.has_g_quadruplex) {
    suggestions.push(
      "G-quadruplex motif detected. These structures are associated with strong target binding in many aptamer families."
    );
  }
  if (features.length < 15 || features.length > 80) {
    suggestions.push(
      "Sequence length is unusual. Most validated aptamers are 20–60 nt."
    );
  }
  if (suggestions.length === 0) {
    suggestions.push(
      "No red flags detected — this sequence's features fall within typical ranges for validated aptamers."
    );
  }
  return suggestions;
}

function downloadCSV(sequence, features, percentiles) {
  const rows = [
    ["metric", "value", "percentile_vs_database"],
    ["sequence", sequence, ""],
    ["type", features.type, ""],
    ["length_nt", features.length, ""],
    ["gc_content", features.gc_content.toFixed(4), percentiles.gc_content?.toFixed(1) ?? ""],
    ["a_freq", features.a_freq.toFixed(4), ""],
    ["t_freq", features.t_freq.toFixed(4), ""],
    ["g_freq", features.g_freq.toFixed(4), ""],
    ["c_freq", features.c_freq.toFixed(4), ""],
    ["purine_ratio", features.purine_ratio.toFixed(4), percentiles.purine_ratio?.toFixed(1) ?? ""],
    ["complexity_shannon_entropy", features.complexity.toFixed(4), percentiles.complexity?.toFixed(1) ?? ""],
    ["has_g_quadruplex", features.has_g_quadruplex, ""],
    ["longest_repeat_nt", features.longest_repeat, ""],
  ];
  const csv = rows.map((r) => r.map((v) => `"${v}"`).join(",")).join("\n");
  const blob = new Blob([csv], { type: "text/csv" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = "aptascope_analysis.csv";
  a.click();
  URL.revokeObjectURL(url);
}

export default function Analyzer({ data }) {
  const [inputText, setInputText] = useState("");
  const [analysis, setAnalysis] = useState(null);
  const [similarResults, setSimilarResults] = useState(null);
  const [isSearching, setIsSearching] = useState(false);
  const [error, setError] = useState(null);
  const workerRef = useRef(null);

  const distributions = useMemo(() => buildDistributions(data), [data]);

  const highAffinityMedians = useMemo(() => {
    const strong = data.filter((r) => r.kd_nM !== null && r.kd_nM < 10);
    const pick = (key) => median(strong.map((r) => r[key]));
    return {
      gc_content: pick("gc_content"),
      complexity: pick("complexity"),
      length: pick("length"),
      purine_ratio: pick("purine_ratio"),
      predicted_mfe: median(strong.filter((r) => r.predicted_mfe !== null).map((r) => r.predicted_mfe)),
      count: strong.length,
    };
  }, [data]);

  useEffect(() => {
    workerRef.current = new Worker(
      new URL("../utils/similarityWorker.js", import.meta.url),
      { type: "module" }
    );
    workerRef.current.onmessage = (e) => {
      setSimilarResults(e.data);
      setIsSearching(false);
    };
    return () => workerRef.current?.terminate();
  }, []);

  function runAnalysis(rawText) {
    const sequence = cleanSequence(rawText);
    if (!sequence) {
      setError("Please paste a valid DNA/RNA sequence (letters A, T, G, C, U only).");
      setAnalysis(null);
      return;
    }
    setError(null);
    const features = computeFeatures(sequence);
    setAnalysis({ sequence, ...features });
    setSimilarResults(null);
    setIsSearching(true);
    workerRef.current.postMessage({ sequence, records: data, topN: 5 });
  }

  function handleExample(example) {
    setInputText(example.sequence);
    runAnalysis(example.sequence);
  }

  const percentiles = analysis
    ? {
        gc_content: percentileRank(analysis.gc_content, distributions.gc_content),
        complexity: percentileRank(analysis.complexity, distributions.complexity),
        length: percentileRank(analysis.length, distributions.length),
        purine_ratio: percentileRank(analysis.purine_ratio, distributions.purine_ratio),
      }
    : {};

  const radarData = analysis
    ? [
        {
          axis: "GC Content",
          user: percentileRank(analysis.gc_content, distributions.gc_content),
          reference: percentileRank(highAffinityMedians.gc_content, distributions.gc_content),
        },
        {
          axis: "Complexity",
          user: percentileRank(analysis.complexity, distributions.complexity),
          reference: percentileRank(highAffinityMedians.complexity, distributions.complexity),
        },
        {
          axis: "Length",
          user: percentileRank(analysis.length, distributions.length),
          reference: percentileRank(highAffinityMedians.length, distributions.length),
        },
        {
          axis: "Purine Ratio",
          user: percentileRank(analysis.purine_ratio, distributions.purine_ratio),
          reference: percentileRank(highAffinityMedians.purine_ratio, distributions.purine_ratio),
        },
      ]
    : [];

  const suggestions = analysis ? generateSuggestions(analysis) : [];

  return (
    <div className="space-y-6">
      <PageHeading>Analyzer</PageHeading>

      <div className="bg-surface border border-border rounded-md p-4 space-y-3">
        <label className="block text-sm font-medium">
          Paste your aptamer sequence (DNA or RNA)
        </label>
        <textarea
          value={inputText}
          onChange={(e) => setInputText(e.target.value)}
          rows={4}
          placeholder="e.g. GGTTGGTGTGGTTGG"
          className="w-full bg-bg border border-border rounded px-3 py-2 font-mono text-sm tracking-[0.05em] focus:outline-none focus:border-accent"
        />
        {inputText.trim() && (
          <div className="text-xs text-textsecondary">
            Detected type:{" "}
            <span className="text-textprimary">
              {inputText.toUpperCase().includes("U") ? "RNA" : "DNA"}
            </span>
          </div>
        )}
        {error && <div className="text-xs text-danger">{error}</div>}
        <div className="flex items-center gap-3">
          <button
            onClick={() => runAnalysis(inputText)}
            className="px-4 py-1.5 bg-accent text-bg font-medium text-sm rounded hover:opacity-90"
          >
            Analyze
          </button>
          <span className="text-xs text-textsecondary">Try an example:</span>
          {EXAMPLES.map((ex) => (
            <button
              key={ex.label}
              onClick={() => handleExample(ex)}
              className="text-xs px-2 py-1 border border-border rounded text-textsecondary hover:text-textprimary hover:border-accent"
            >
              {ex.label}
            </button>
          ))}
        </div>
      </div>

      {analysis && (
        <>
          {/* Panel A: Computed Features */}
          <div className="bg-surface border border-border rounded-md p-4">
            <div className="flex items-center justify-between mb-3">
              <h2 className="text-xs uppercase tracking-wider font-semibold">Computed Features</h2>
              <button
                onClick={() => downloadCSV(analysis.sequence, analysis, percentiles)}
                className="text-xs px-2 py-1 border border-border rounded text-textsecondary hover:text-textprimary hover:border-accent"
              >
                Export CSV
              </button>
            </div>

            <div className="mb-3">
              <ColoredSequence sequence={analysis.sequence} className="text-sm" />
              <CopyButton text={analysis.sequence} className="ml-2" />
            </div>

            <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
              <FeatureCard label="Length" value={`${analysis.length} nt`} />
              <FeatureCard
                label="GC Content"
                value={`${(analysis.gc_content * 100).toFixed(1)}%`}
                percentile={percentiles.gc_content}
              />
              <FeatureCard
                label="Complexity"
                value={analysis.complexity.toFixed(2)}
                percentile={percentiles.complexity}
              />
              <FeatureCard
                label="G-Quadruplex"
                value={analysis.has_g_quadruplex ? "Yes" : "No"}
              />
              <FeatureCard label="Longest Repeat" value={`${analysis.longest_repeat} nt`} />
              <FeatureCard
                label="Purine Ratio"
                value={`${(analysis.purine_ratio * 100).toFixed(1)}%`}
                percentile={percentiles.purine_ratio}
              />
            </div>
          </div>

          <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
            {/* Panel B: Radar Chart */}
            <div className="bg-surface border border-border rounded-md p-4">
              <h2 className="text-xs uppercase tracking-wider font-semibold mb-1">
                Comparison to High-Affinity Binders
              </h2>
              <p className="text-[13px] text-textsecondary mb-2">
                Your sequence (rose) vs. median of known Kd &lt; 10 nM binders
                (green, n={highAffinityMedians.count})
              </p>
              <div style={{ width: "100%", height: 300 }}>
                <ResponsiveContainer>
                  <RadarChart data={radarData}>
                    <PolarGrid stroke="#2F343C" />
                    <PolarAngleAxis dataKey="axis" tick={{ fontSize: 11, fill: "#8F99A8" }} />
                    <PolarRadiusAxis domain={[0, 100]} tick={{ fontSize: 10, fill: "#8F99A8" }} />
                    <Radar
                      name="Your sequence"
                      dataKey="user"
                      stroke="#E0576B"
                      fill="#E0576B"
                      fillOpacity={0.35}
                    />
                    <Radar
                      name="High-affinity median"
                      dataKey="reference"
                      stroke="#238551"
                      fill="#238551"
                      fillOpacity={0.2}
                    />
                    <Legend wrapperStyle={{ fontSize: 11 }} />
                    <Tooltip
                      formatter={(v) => (v === null ? "N/A" : `${v.toFixed(0)}th pct`)}
                    />
                  </RadarChart>
                </ResponsiveContainer>
              </div>
              <p className="text-[13px] text-textsecondary">
                Axes show each feature's percentile rank within the full database.
              </p>
            </div>

            {/* Panel D: Design Suggestions */}
            <div className="bg-surface border border-border rounded-md p-4">
              <h2 className="text-xs uppercase tracking-wider font-semibold mb-3">Design Suggestions</h2>
              <ul className="space-y-2">
                {suggestions.map((s, i) => (
                  <li key={i} className="text-sm text-textsecondary flex gap-2">
                    <span className="text-accent not-italic font-sans">•</span>
                    <span>{s}</span>
                  </li>
                ))}
              </ul>
            </div>
          </div>

          {/* Panel C: Most Similar Known Aptamers */}
          <div className="bg-surface border border-border rounded-md p-4">
            <h2 className="text-xs uppercase tracking-wider font-semibold mb-3">Most Similar Known Aptamers</h2>
            {isSearching ? (
              <div className="text-sm text-textsecondary py-6 text-center">
                Searching {data.length.toLocaleString()} sequences…
              </div>
            ) : (
              <div className="overflow-x-auto">
                <table className="w-full text-sm">
                  <thead>
                    <tr className="text-left text-textsecondary border-b border-border uppercase tracking-wider text-[11px]">
                      <th className="px-2 py-1.5">Target</th>
                      <th className="px-2 py-1.5">Sequence</th>
                      <th className="px-2 py-1.5">Kd</th>
                      <th className="px-2 py-1.5">Similarity</th>
                      <th className="px-2 py-1.5">Source</th>
                    </tr>
                  </thead>
                  <tbody>
                    {similarResults?.map((r) => (
                      <tr key={r.id} className="border-b border-border/60">
                        <td className="px-2 py-1.5">{r.target_name}</td>
                        <td className="px-2 py-1.5">
                          <div className="flex items-center gap-2">
                            <ColoredSequence
                              sequence={
                                r.sequence.length > 40
                                  ? r.sequence.slice(0, 40) + "…"
                                  : r.sequence
                              }
                              className="whitespace-nowrap"
                            />
                            <CopyButton text={r.sequence} />
                          </div>
                        </td>
                        <td className={`px-2 py-1.5 font-mono ${kdColorClass(r.kd_nM)}`}>
                          {formatKd(r.kd_nM)}
                        </td>
                        <td className="px-2 py-1.5">{r.similarity.toFixed(1)}%</td>
                        <td className="px-2 py-1.5 text-textsecondary">{r.source}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            )}
          </div>
        </>
      )}
    </div>
  );
}

function FeatureCard({ label, value, percentile }) {
  return (
    <div className="bg-bg border border-border rounded p-3">
      <div className="text-xs uppercase tracking-wider text-textsecondary mb-1">{label}</div>
      <div className="text-lg font-mono font-medium">{value}</div>
      {percentile !== undefined && percentile !== null && (
        <>
          <div className="text-xs text-textsecondary mt-1 mb-1">
            {ordinalSuffix(percentile)} percentile
          </div>
          <PercentileGauge percentile={percentile} />
        </>
      )}
    </div>
  );
}
