import { useEffect, useRef, useState } from "react";
// See ProteinStructureViewer.jsx for why this is a bare side-effect import.
import "3dmol";

// Renders a real, experimentally-solved RCSB structure — as opposed to
// ProteinStructureViewer's AlphaFold prediction, this can show the
// aptamer's own 3D fold (AlphaFold has no model for nucleic acids at all).
export default function PdbStructureViewer({ entry }) {
  const containerRef = useRef(null);
  const viewerRef = useRef(null);
  const [pdbData, setPdbData] = useState(null);
  const [status, setStatus] = useState("loading"); // loading | found | error

  useEffect(() => {
    let cancelled = false;
    setStatus("loading");
    setPdbData(null);

    fetch(`https://files.rcsb.org/download/${entry.pdb_id}.pdb`)
      .then((res) => {
        if (!res.ok) throw new Error(res.statusText);
        return res.text();
      })
      .then((data) => {
        if (cancelled) return;
        setPdbData(data);
        setStatus("found");
      })
      .catch(() => {
        if (!cancelled) setStatus("error");
      });

    return () => {
      cancelled = true;
    };
  }, [entry.pdb_id]);

  useEffect(() => {
    if (status !== "found" || !pdbData || !containerRef.current) return;

    containerRef.current.innerHTML = "";
    const viewer = window.$3Dmol.createViewer(containerRef.current, {
      backgroundColor: "#111418",
    });
    viewerRef.current = viewer;
    viewer.addModel(pdbData, "pdb");
    // cartoon covers the polymer chain (protein ribbon or nucleic acid
    // tube — 3Dmol auto-detects which); stick covers bound ligands/ions
    // (HETATM groups), several of these entries are aptamer-ligand
    // complexes (e.g. malachite green, theophylline, DFHBI) that cartoon
    // alone would render as nothing at all.
    viewer.setStyle({}, { cartoon: { color: "spectrum" } });
    viewer.setStyle({ hetflag: true }, { stick: {} });
    viewer.zoomTo();
    viewer.render();
    viewer.spin("y", 0.4);

    return () => {
      viewerRef.current?.clear();
      viewerRef.current = null;
    };
  }, [status, pdbData]);

  if (status === "loading") {
    return (
      <div className="text-sm text-textsecondary py-8 text-center">
        Loading {entry.pdb_id} from RCSB PDB…
      </div>
    );
  }

  if (status === "error") {
    return (
      <div className="text-sm text-danger py-8 text-center">
        Couldn't load {entry.pdb_id} from RCSB — try again later.
      </div>
    );
  }

  return (
    <div>
      <div ref={containerRef} style={{ width: "100%", height: 320, position: "relative" }} />
      <div className="flex items-center justify-between mt-2 text-xs text-textsecondary">
        <span>
          {entry.title}
          {" · "}
          <span className="font-mono">{entry.method}</span>
          {entry.resolution && (
            <>
              {" · "}
              <span className="font-mono">{entry.resolution} Å</span>
            </>
          )}
        </span>
        <a
          href={entry.viewer_url || `https://www.rcsb.org/structure/${entry.pdb_id}`}
          target="_blank"
          rel="noreferrer"
          className="text-accent hover:underline shrink-0 ml-3"
        >
          View on RCSB PDB ↗
        </a>
      </div>
    </div>
  );
}
