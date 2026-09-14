import { Fragment, useMemo, useState } from "react";
import ColoredSequence, { CopyButton } from "./ColoredSequence.jsx";
import SequenceDetail from "./SequenceDetail.jsx";
import PageHeading from "./PageHeading.jsx";
import { formatKd, kdColorClass } from "../utils/sequence.js";
import pdbEnrichment from "../data/pdb_enrichment.json";

const TARGET_TYPE_OPTIONS = ["All", "Protein", "Small Molecule", "Cell", "Nucleic Acid", "Microorganism", "Other"];

const KD_LOG_MIN = -2; // 0.01 nM
const KD_LOG_MAX = 5; // 100,000 nM
const LENGTH_MIN = 10;
const LENGTH_MAX = 200;

const PAGE_SIZE = 50;

function truncateSeq(seq, n = 40) {
  return seq.length > n ? seq.slice(0, n) + "…" : seq;
}

export default function Explorer({ data }) {
  const [search, setSearch] = useState("");
  const [typeFilter, setTypeFilter] = useState("All");
  const [targetTypeFilter, setTargetTypeFilter] = useState("All");
  const [has3dOnly, setHas3dOnly] = useState(false);
  const [kdLogRange, setKdLogRange] = useState([KD_LOG_MIN, KD_LOG_MAX]);
  const [lengthRange, setLengthRange] = useState([LENGTH_MIN, LENGTH_MAX]);
  const [sortKey, setSortKey] = useState("kd_nM");
  const [sortDir, setSortDir] = useState("asc");
  const [expandedId, setExpandedId] = useState(null);
  const [page, setPage] = useState(1);

  const kdMin = 10 ** kdLogRange[0];
  const kdMax = 10 ** kdLogRange[1];

  // Target names known to have a solved 3D structure — approximated by
  // substring-matching each distinct target name against RCSB PDB entry
  // titles, since none of the source databases carry a per-record PDB ID.
  const structureTargets = useMemo(() => {
    const titles = Object.values(pdbEnrichment).map((p) => p.title.toLowerCase());
    const set = new Set();
    for (const name of new Set(data.map((r) => r.target_name))) {
      const lower = name.toLowerCase();
      if (lower.length >= 4 && titles.some((t) => t.includes(lower))) {
        set.add(name);
      }
    }
    return set;
  }, [data]);

  const filtered = useMemo(() => {
    const q = search.trim().toLowerCase();
    let rows = data.filter((r) => {
      if (typeFilter !== "All" && r.aptamer_type !== typeFilter) return false;
      if (targetTypeFilter !== "All" && r.target_type !== targetTypeFilter) return false;
      if (has3dOnly && !structureTargets.has(r.target_name)) return false;
      if (r.length < lengthRange[0] || r.length > lengthRange[1]) return false;
      if (r.kd_nM !== null && (r.kd_nM < kdMin || r.kd_nM > kdMax)) return false;
      if (q) {
        const hay = `${r.target_name} ${r.sequence}`.toLowerCase();
        if (!hay.includes(q)) return false;
      }
      return true;
    });

    rows = [...rows].sort((a, b) => {
      let av = a[sortKey];
      let bv = b[sortKey];
      if (av === null || av === undefined) av = sortDir === "asc" ? Infinity : -Infinity;
      if (bv === null || bv === undefined) bv = sortDir === "asc" ? Infinity : -Infinity;
      if (typeof av === "string") {
        return sortDir === "asc" ? av.localeCompare(bv) : bv.localeCompare(av);
      }
      return sortDir === "asc" ? av - bv : bv - av;
    });

    return rows;
  }, [data, search, typeFilter, targetTypeFilter, has3dOnly, structureTargets, kdLogRange, lengthRange, sortKey, sortDir]);

  const totalPages = Math.max(1, Math.ceil(filtered.length / PAGE_SIZE));
  const pageRows = filtered.slice((page - 1) * PAGE_SIZE, page * PAGE_SIZE);

  function handleSort(key) {
    if (sortKey === key) {
      setSortDir(sortDir === "asc" ? "desc" : "asc");
    } else {
      setSortKey(key);
      setSortDir("asc");
    }
    setPage(1);
  }

  function resetPage(setter) {
    return (val) => {
      setter(val);
      setPage(1);
    };
  }

  const columns = [
    { key: "target_name", label: "Target" },
    { key: "aptamer_type", label: "Type" },
    { key: "sequence", label: "Sequence", sortable: false },
    { key: "length", label: "Length" },
    { key: "kd_nM", label: "Kd (nM)" },
    { key: "gc_content", label: "GC%" },
  ];

  return (
    <div className="space-y-4">
      <PageHeading>Explorer</PageHeading>

      <div className="bg-surface border border-border rounded-md p-4 grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-4">
        <div className="lg:col-span-2">
          <label className="block text-xs uppercase tracking-wider text-textsecondary mb-1">Search (target or sequence)</label>
          <input
            type="text"
            value={search}
            onChange={(e) => resetPage(setSearch)(e.target.value)}
            placeholder="e.g. thrombin, GGTTGG…"
            className="w-full bg-bg border border-border rounded px-3 py-1.5 text-sm focus:outline-none focus:border-accent"
          />
        </div>

        <div>
          <label className="block text-xs uppercase tracking-wider text-textsecondary mb-1">Type</label>
          <div className="flex gap-1">
            {["All", "DNA", "RNA"].map((t) => (
              <button
                key={t}
                onClick={() => resetPage(setTypeFilter)(t)}
                className={`flex-1 text-sm px-2 py-1.5 rounded-full border ${
                  typeFilter === t
                    ? "border-accent text-accent"
                    : "border-border text-textsecondary hover:text-textprimary"
                }`}
              >
                {t}
              </button>
            ))}
          </div>
        </div>

        <div>
          <label className="block text-xs uppercase tracking-wider text-textsecondary mb-1">Target Type</label>
          <select
            value={targetTypeFilter}
            onChange={(e) => resetPage(setTargetTypeFilter)(e.target.value)}
            className="w-full bg-bg border border-border rounded px-3 py-1.5 text-sm focus:outline-none focus:border-accent"
          >
            {TARGET_TYPE_OPTIONS.map((label) => (
              <option key={label} value={label}>
                {label}
              </option>
            ))}
          </select>
        </div>

        <div>
          <label className="block text-xs uppercase tracking-wider text-textsecondary mb-1">Structure</label>
          <button
            onClick={() => resetPage(setHas3dOnly)(!has3dOnly)}
            className={`w-full text-sm px-2 py-1.5 rounded-full border ${
              has3dOnly
                ? "border-accent text-accent"
                : "border-border text-textsecondary hover:text-textprimary"
            }`}
          >
            Has 3D Structure
          </button>
        </div>

        <div className="lg:col-span-2">
          <label className="block text-xs text-textsecondary mb-1">
            <span className="uppercase tracking-wider">Kd range</span>: {formatKd(kdMin)} – {formatKd(kdMax)}
          </label>
          <div className="flex gap-2 items-center">
            <input
              type="range"
              min={KD_LOG_MIN}
              max={KD_LOG_MAX}
              step={0.1}
              value={kdLogRange[0]}
              onChange={(e) =>
                resetPage(setKdLogRange)([
                  Math.min(Number(e.target.value), kdLogRange[1]),
                  kdLogRange[1],
                ])
              }
              className="flex-1 accent-accent"
            />
            <input
              type="range"
              min={KD_LOG_MIN}
              max={KD_LOG_MAX}
              step={0.1}
              value={kdLogRange[1]}
              onChange={(e) =>
                resetPage(setKdLogRange)([
                  kdLogRange[0],
                  Math.max(Number(e.target.value), kdLogRange[0]),
                ])
              }
              className="flex-1 accent-accent"
            />
          </div>
        </div>

        <div className="lg:col-span-2">
          <label className="block text-xs text-textsecondary mb-1">
            <span className="uppercase tracking-wider">Length range</span>: {lengthRange[0]} – {lengthRange[1]} nt
          </label>
          <div className="flex gap-2 items-center">
            <input
              type="range"
              min={LENGTH_MIN}
              max={LENGTH_MAX}
              value={lengthRange[0]}
              onChange={(e) =>
                resetPage(setLengthRange)([
                  Math.min(Number(e.target.value), lengthRange[1]),
                  lengthRange[1],
                ])
              }
              className="flex-1 accent-accent"
            />
            <input
              type="range"
              min={LENGTH_MIN}
              max={LENGTH_MAX}
              value={lengthRange[1]}
              onChange={(e) =>
                resetPage(setLengthRange)([
                  lengthRange[0],
                  Math.max(Number(e.target.value), lengthRange[0]),
                ])
              }
              className="flex-1 accent-accent"
            />
          </div>
        </div>
      </div>

      <div className="text-xs text-textsecondary">
        {filtered.length.toLocaleString()} of {data.length.toLocaleString()} aptamers
      </div>

      <div className="bg-surface border border-border rounded-md overflow-x-auto">
        <table className="w-full text-sm">
          <thead>
            <tr className="border-b border-border text-left text-textsecondary bg-bg/40">
              {columns.map((col) => (
                <th
                  key={col.key}
                  onClick={() => col.sortable !== false && handleSort(col.key)}
                  className={`px-3 py-2 font-medium uppercase tracking-wider text-[11px] whitespace-nowrap ${
                    col.sortable !== false ? "cursor-pointer hover:text-textprimary" : ""
                  }`}
                >
                  {col.label}
                  {sortKey === col.key && (sortDir === "asc" ? " ↑" : " ↓")}
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            {pageRows.map((r) => (
              <Fragment key={r.id}>
                <tr
                  onClick={() => setExpandedId(expandedId === r.id ? null : r.id)}
                  className="border-b border-border/60 hover:bg-bg cursor-pointer"
                >
                  <td className="px-3 py-2">{r.target_name}</td>
                  <td className="px-3 py-2">
                    <span
                      className={`text-xs px-2 py-0.5 rounded-full border ${
                        r.aptamer_type === "DNA"
                          ? "border-dna text-dna"
                          : "border-rna text-rna"
                      }`}
                    >
                      {r.aptamer_type}
                    </span>
                  </td>
                  <td className="px-3 py-2">
                    <div className="flex items-center gap-2">
                      <ColoredSequence sequence={truncateSeq(r.sequence)} className="whitespace-nowrap" />
                      <CopyButton text={r.sequence} />
                    </div>
                  </td>
                  <td className="px-3 py-2 font-mono">{r.length}</td>
                  <td className={`px-3 py-2 font-mono ${kdColorClass(r.kd_nM)}`}>
                    {formatKd(r.kd_nM)}
                  </td>
                  <td className="px-3 py-2 font-mono">{(r.gc_content * 100).toFixed(1)}%</td>
                </tr>
                {expandedId === r.id && (
                  <tr>
                    <td colSpan={columns.length} className="p-3 bg-bg/40">
                      <SequenceDetail record={r} />
                    </td>
                  </tr>
                )}
              </Fragment>
            ))}
          </tbody>
        </table>
      </div>

      {totalPages > 1 && (
        <div className="flex items-center justify-center gap-3 text-sm">
          <button
            disabled={page <= 1}
            onClick={() => setPage((p) => Math.max(1, p - 1))}
            className="px-3 py-1 border border-border rounded disabled:opacity-40 hover:border-accent"
          >
            Prev
          </button>
          <span className="text-textsecondary">
            Page {page} of {totalPages}
          </span>
          <button
            disabled={page >= totalPages}
            onClick={() => setPage((p) => Math.min(totalPages, p + 1))}
            className="px-3 py-1 border border-border rounded disabled:opacity-40 hover:border-accent"
          >
            Next
          </button>
        </div>
      )}
    </div>
  );
}
