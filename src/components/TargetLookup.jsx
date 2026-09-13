import { useMemo, useState } from "react";
import {
  ResponsiveContainer,
  ScatterChart,
  Scatter,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
} from "recharts";
import ColoredSequence, { CopyButton } from "./ColoredSequence.jsx";
import PageHeading from "./PageHeading.jsx";
import StructureViewer from "./StructureViewer.jsx";
import ProteinStructureViewer from "./ProteinStructureViewer.jsx";
import { formatKd, kdColorClass } from "../utils/sequence.js";
import { quantile } from "../utils/stats.js";
import { resolveUniProtId } from "../utils/alphafold.js";

const PAGE_SIZE = 50;

function TooltipBox({ active, payload }) {
  if (!active || !payload || !payload.length) return null;
  const p = payload[0].payload;
  return (
    <div className="bg-bg border border-border rounded px-2 py-1.5 text-xs">
      <div>GC: {(p.x * 100).toFixed(1)}%</div>
      <div>Kd: {formatKd(p.y)}</div>
    </div>
  );
}

export default function TargetLookup({ data }) {
  const targets = useMemo(() => {
    const set = new Set(data.map((r) => r.target_name));
    return [...set].sort((a, b) => a.localeCompare(b));
  }, [data]);

  const [selected, setSelected] = useState("");
  const [query, setQuery] = useState("");
  const [sortKey, setSortKey] = useState("kd_nM");
  const [sortDir, setSortDir] = useState("asc");
  const [page, setPage] = useState(1);

  const filteredTargets = useMemo(() => {
    const q = query.trim().toLowerCase();
    if (!q) return targets;
    return targets.filter((t) => t.toLowerCase().includes(q));
  }, [targets, query]);

  const records = useMemo(
    () => (selected ? data.filter((r) => r.target_name === selected) : []),
    [data, selected]
  );

  const sortedRecords = useMemo(() => {
    return [...records].sort((a, b) => {
      let av = a[sortKey];
      let bv = b[sortKey];
      if (av === null || av === undefined) av = sortDir === "asc" ? Infinity : -Infinity;
      if (bv === null || bv === undefined) bv = sortDir === "asc" ? Infinity : -Infinity;
      return sortDir === "asc" ? av - bv : bv - av;
    });
  }, [records, sortKey, sortDir]);

  const totalPages = Math.max(1, Math.ceil(sortedRecords.length / PAGE_SIZE));
  const pageRecords = sortedRecords.slice((page - 1) * PAGE_SIZE, page * PAGE_SIZE);

  function handleSort(key) {
    if (sortKey === key) setSortDir(sortDir === "asc" ? "desc" : "asc");
    else {
      setSortKey(key);
      setSortDir("asc");
    }
    setPage(1);
  }

  function handleSelectTarget(t) {
    setSelected(t);
    setPage(1);
  }

  const profile = useMemo(() => {
    if (records.length === 0) return null;
    const gcs = records.map((r) => r.gc_content).sort((a, b) => a - b);
    const lengths = records.map((r) => r.length).sort((a, b) => a - b);
    const kds = records.map((r) => r.kd_nM).filter((k) => k !== null);
    const bestKd = kds.length ? Math.min(...kds) : null;

    const gqCount = records.filter((r) => r.has_g_quadruplex).length;
    const gqPct = (gqCount / records.length) * 100;

    const avgFreq = (key) => records.reduce((s, r) => s + r[key], 0) / records.length;
    const overallAvg = { A: 0.25, T: 0.25, G: 0.25, C: 0.25 }; // uniform baseline
    const freqs = { A: avgFreq("a_freq"), T: avgFreq("t_freq"), G: avgFreq("g_freq"), C: avgFreq("c_freq") };
    let biasedBase = null;
    let maxDev = 0.05; // minimum deviation to call out a bias
    for (const base of ["A", "T", "G", "C"]) {
      const dev = freqs[base] - overallAvg[base];
      if (dev > maxDev) {
        maxDev = dev;
        biasedBase = base;
      }
    }

    return {
      bestKd,
      gcRange: [quantile(gcs, 0.25), quantile(gcs, 0.75)],
      lengthRange: [quantile(lengths, 0.25), quantile(lengths, 0.75)],
      gqPct,
      biasedBase,
      scatter: records.filter((r) => r.kd_nM !== null).map((r) => ({ x: r.gc_content, y: r.kd_nM })),
    };
  }, [records]);

  const topBinders = useMemo(() => {
    return [...records]
      .filter((r) => r.kd_nM !== null && r.mfe_structure)
      .sort((a, b) => a.kd_nM - b.kd_nM)
      .slice(0, 3);
  }, [records]);

  const uniprotId = useMemo(() => resolveUniProtId(records), [records]);

  return (
    <div className="space-y-4">
      <PageHeading>Target Lookup</PageHeading>

      <div className="bg-surface border border-border rounded-md p-4 space-y-2">
        <label className="block text-xs uppercase tracking-wider text-textsecondary">Select a target</label>
        <input
          type="text"
          value={query}
          onChange={(e) => setQuery(e.target.value)}
          placeholder="Search targets…"
          className="w-full bg-bg border border-border rounded px-3 py-1.5 text-sm focus:outline-none focus:border-accent mb-2"
        />
        <select
          value={selected}
          onChange={(e) => handleSelectTarget(e.target.value)}
          size={Math.min(8, Math.max(4, filteredTargets.length))}
          className="w-full bg-bg border border-border rounded px-3 py-1.5 text-sm focus:outline-none focus:border-accent"
        >
          {filteredTargets.map((t) => (
            <option key={t} value={t}>
              {t}
            </option>
          ))}
        </select>
      </div>

      {selected && profile && (
        <>
          <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
            <StatCard label="Known Aptamers" value={records.length} />
            <StatCard
              label="Best Kd"
              value={formatKd(profile.bestKd)}
              valueClass={kdColorClass(profile.bestKd)}
            />
            <StatCard
              label="Optimal GC Range"
              value={`${(profile.gcRange[0] * 100).toFixed(0)}–${(profile.gcRange[1] * 100).toFixed(0)}%`}
            />
            <StatCard
              label="Optimal Length Range"
              value={`${profile.lengthRange[0].toFixed(0)}–${profile.lengthRange[1].toFixed(0)} nt`}
            />
          </div>

          <div className="bg-surface border border-border rounded-md p-4">
            <h2 className="text-xs uppercase tracking-wider font-semibold mb-3">Design Profile — {selected}</h2>
            <div className="grid grid-cols-1 md:grid-cols-2 gap-4 text-sm text-textprimary">
              <div>
                <div className="text-xs uppercase tracking-wider text-textsecondary mb-0.5 font-sans">
                  Common Structural Motifs
                </div>
                <div>
                  {profile.gqPct > 15
                    ? `G-quadruplex motif present in ${profile.gqPct.toFixed(0)}% of known binders`
                    : "No strongly recurring motif detected"}
                </div>
              </div>
              <div>
                <div className="text-xs uppercase tracking-wider text-textsecondary mb-0.5 font-sans">
                  Common Nucleotide Bias
                </div>
                <div>
                  {profile.biasedBase
                    ? `${profile.biasedBase}-rich sequences are over-represented`
                    : "No strong nucleotide bias detected"}
                </div>
              </div>
            </div>
          </div>

          {uniprotId && (
            <div className="bg-surface border border-border rounded-md p-4">
              <h2 className="text-xs uppercase tracking-wider font-semibold mb-3">
                Target 3D Structure (AlphaFold)
              </h2>
              <p className="text-[13px] text-textsecondary mb-3">
                This is the predicted structure of {selected} itself, the protein this
                aptamer binds — not the aptamer's own fold, which AlphaFold can't
                predict (it doesn't model DNA/RNA).
              </p>
              <ProteinStructureViewer uniprotId={uniprotId} targetName={selected} />
            </div>
          )}

          {topBinders.length > 0 && (
            <div className="bg-surface border border-border rounded-md p-4">
              <h2 className="text-xs uppercase tracking-wider font-semibold mb-3">
                Top Binder Structures
              </h2>
              <div
                className={`grid grid-cols-1 gap-4 ${
                  { 1: "md:grid-cols-1", 2: "md:grid-cols-2", 3: "md:grid-cols-3" }[topBinders.length]
                }`}
              >
                {topBinders.map((r) => (
                  <div key={r.id}>
                    <StructureViewer
                      sequence={r.sequence}
                      dotBracket={r.mfe_structure}
                      label={`Kd = ${formatKd(r.kd_nM)} · ${r.length} nt`}
                      height={220}
                    />
                  </div>
                ))}
              </div>
            </div>
          )}

          <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
            <div className="bg-surface border border-border rounded-md p-4">
              <h2 className="text-xs uppercase tracking-wider font-semibold mb-3">GC Content vs. Kd</h2>
              <div style={{ width: "100%", height: 280 }}>
                <ResponsiveContainer>
                  <ScatterChart margin={{ top: 10, right: 10, bottom: 10, left: 0 }}>
                    <CartesianGrid stroke="#2F343C" strokeDasharray="3 3" />
                    <XAxis
                      type="number"
                      dataKey="x"
                      domain={[0, 1]}
                      tick={{ fontSize: 11, fill: "#8F99A8" }}
                      label={{ value: "GC Content", position: "insideBottom", offset: -5, style: { fontSize: 11, fill: "#8F99A8" } }}
                    />
                    <YAxis
                      type="number"
                      dataKey="y"
                      scale="log"
                      domain={["auto", "auto"]}
                      tick={{ fontSize: 11, fill: "#8F99A8" }}
                      label={{ value: "Kd (nM, log)", angle: -90, position: "insideLeft", style: { fontSize: 11, fill: "#8F99A8" } }}
                    />
                    <Tooltip content={<TooltipBox />} />
                    <Scatter data={profile.scatter} fill="#E0576B" />
                  </ScatterChart>
                </ResponsiveContainer>
              </div>
            </div>

            <div className="bg-surface border border-border rounded-md p-4 overflow-x-auto">
              <h2 className="text-xs uppercase tracking-wider font-semibold mb-3">
                All Aptamers for {selected} ({records.length.toLocaleString()})
              </h2>
              <table className="w-full text-sm">
                <thead>
                  <tr className="text-left text-textsecondary border-b border-border uppercase tracking-wider text-[11px]">
                    <th
                      className="px-2 py-1.5 cursor-pointer hover:text-textprimary"
                      onClick={() => handleSort("sequence")}
                    >
                      Sequence
                    </th>
                    <th
                      className="px-2 py-1.5 cursor-pointer hover:text-textprimary"
                      onClick={() => handleSort("kd_nM")}
                    >
                      Kd
                    </th>
                    <th
                      className="px-2 py-1.5 cursor-pointer hover:text-textprimary"
                      onClick={() => handleSort("length")}
                    >
                      Length
                    </th>
                  </tr>
                </thead>
                <tbody>
                  {pageRecords.map((r) => (
                    <tr key={r.id} className="border-b border-border/60">
                      <td className="px-2 py-1.5">
                        <div className="flex items-center gap-2">
                          <ColoredSequence
                            sequence={r.sequence.length > 30 ? r.sequence.slice(0, 30) + "…" : r.sequence}
                            className="whitespace-nowrap"
                          />
                          <CopyButton text={r.sequence} />
                        </div>
                      </td>
                      <td className={`px-2 py-1.5 font-mono ${kdColorClass(r.kd_nM)}`}>
                        {formatKd(r.kd_nM)}
                      </td>
                      <td className="px-2 py-1.5 font-mono">{r.length}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
              {totalPages > 1 && (
                <div className="flex items-center justify-center gap-3 text-sm mt-3 pt-3 border-t border-border">
                  <button
                    disabled={page <= 1}
                    onClick={() => setPage((p) => Math.max(1, p - 1))}
                    className="px-3 py-1 border border-border disabled:opacity-40 hover:border-accent"
                  >
                    Prev
                  </button>
                  <span className="text-textsecondary font-mono">
                    Page {page} of {totalPages}
                  </span>
                  <button
                    disabled={page >= totalPages}
                    onClick={() => setPage((p) => Math.min(totalPages, p + 1))}
                    className="px-3 py-1 border border-border disabled:opacity-40 hover:border-accent"
                  >
                    Next
                  </button>
                </div>
              )}
            </div>
          </div>
        </>
      )}
    </div>
  );
}

function StatCard({ label, value, valueClass = "" }) {
  return (
    <div className="bg-surface border border-border rounded-md p-4">
      <div className="text-xs uppercase tracking-wider text-textsecondary mb-1">{label}</div>
      <div className={`text-lg font-mono font-medium ${valueClass}`}>{value}</div>
    </div>
  );
}
