import { useState } from "react";
import PdbStructureViewer from "./PdbStructureViewer.jsx";

// A target can match more than one solved structure (e.g. several crystal
// forms of the same aptamer) — offer a picker when there's more than one.
export default function PdbStructurePanel({ matches }) {
  const [index, setIndex] = useState(0);
  const entry = matches[Math.min(index, matches.length - 1)];

  return (
    <div>
      {matches.length > 1 && (
        <div className="flex gap-1.5 mb-3 flex-wrap">
          {matches.map((m, i) => (
            <button
              key={m.pdb_id}
              onClick={() => setIndex(i)}
              className={`text-xs font-mono px-2 py-1 rounded-full border ${
                i === index
                  ? "border-accent text-accent"
                  : "border-border text-textsecondary hover:text-textprimary"
              }`}
            >
              {m.pdb_id}
            </button>
          ))}
        </div>
      )}
      <PdbStructureViewer key={entry.pdb_id} entry={entry} />
    </div>
  );
}
