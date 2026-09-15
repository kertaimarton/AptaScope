import { useMemo } from "react";
import ColoredSequence, { CopyButton } from "./ColoredSequence.jsx";
import StructureViewer from "./StructureViewer.jsx";
import PdbStructurePanel from "./PdbStructurePanel.jsx";
import { formatKd, kdColorClass } from "../utils/sequence.js";
import { findPdbStructures } from "../utils/pdb.js";

const FEATURE_ROWS = [
  ["GC Content", (r) => `${(r.gc_content * 100).toFixed(1)}%`],
  ["A / T(U) / G / C", (r) =>
    `${(r.a_freq * 100).toFixed(0)}% / ${(r.t_freq * 100).toFixed(0)}% / ${(r.g_freq * 100).toFixed(0)}% / ${(r.c_freq * 100).toFixed(0)}%`],
  ["Purine Ratio", (r) => `${(r.purine_ratio * 100).toFixed(1)}%`],
  ["Complexity (Shannon entropy)", (r) => r.complexity.toFixed(3)],
  ["G-Quadruplex Motif", (r) => (r.has_g_quadruplex ? "Yes" : "No")],
  ["Longest Homopolymer Run", (r) => `${r.longest_repeat} nt`],
  ["Predicted MFE", (r) => (r.predicted_mfe !== null ? `${r.predicted_mfe} kcal/mol` : "—")],
  ["Predicted Stems / Loops", (r) =>
    r.num_stems !== null ? `${r.num_stems} / ${r.num_loops}` : "—"],
  ["SELEX Method", (r) => r.selex_method || "Not reported"],
  ["Year", (r) => r.year || "Unknown"],
  ["Source", (r) => r.source],
];

export default function SequenceDetail({ record }) {
  const pdbMatches = useMemo(() => findPdbStructures(record.target_name), [record.target_name]);

  return (
    <div className="bg-bg border border-border rounded-md p-4 space-y-4">
      <div className="flex items-start justify-between gap-4">
        <div className="min-w-0">
          <div className="text-xs uppercase tracking-wider text-textsecondary mb-1">
            Full sequence ({record.length} nt)
          </div>
          <ColoredSequence sequence={record.sequence} className="text-sm" />
        </div>
        <CopyButton text={record.sequence} className="shrink-0" />
      </div>

      <div className="grid grid-cols-2 md:grid-cols-3 gap-3 text-sm">
        {FEATURE_ROWS.map(([label, fn]) => (
          <div key={label}>
            <div className="text-xs uppercase tracking-wider text-textsecondary">{label}</div>
            <div>{fn(record)}</div>
          </div>
        ))}
        <div>
          <div className="text-xs uppercase tracking-wider text-textsecondary">Kd</div>
          <div className={kdColorClass(record.kd_nM)}>{formatKd(record.kd_nM)}</div>
        </div>
      </div>

      {(record.cancer_type || record.bacterial_species) && (
        <div className="flex gap-6 text-sm">
          {record.cancer_type && (
            <div>
              <div className="text-xs uppercase tracking-wider text-textsecondary">Cancer Type</div>
              <div>{record.cancer_type}</div>
            </div>
          )}
          {record.bacterial_species && (
            <div>
              <div className="text-xs uppercase tracking-wider text-textsecondary">Bacterial Species</div>
              <div>{record.bacterial_species}</div>
            </div>
          )}
        </div>
      )}

      {record.doi && (
        <a
          href={`https://doi.org/${record.doi}`}
          target="_blank"
          rel="noreferrer"
          className="inline-block font-mono text-sm text-accent hover:underline"
        >
          View publication (DOI: {record.doi}) ↗
        </a>
      )}

      {pdbMatches.length > 0 && (
        <div className="border-t border-border pt-4">
          <div className="text-xs uppercase tracking-wider text-textsecondary mb-2">
            Solved 3D Structure (RCSB PDB)
          </div>
          <PdbStructurePanel matches={pdbMatches} />
        </div>
      )}

      {record.mfe_structure && (
        <div className="border-t border-border pt-4">
          <StructureViewer
            sequence={record.sequence}
            dotBracket={record.mfe_structure}
            label="Secondary Structure"
            height={260}
          />
        </div>
      )}
    </div>
  );
}
