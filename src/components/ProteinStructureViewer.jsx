import { useEffect, useRef, useState } from "react";
// 3dmol ships as a UMD build with no real ESM/CJS export shape (its default
// export, once Vite's CJS interop unwraps it, ends up wrapped in an extra
// `.default` layer that varies by import form). The one thing that's
// reliably true in a browser is that the UMD wrapper attaches itself to
// `window.$3Dmol` as a side effect — so import it purely for that effect
// and read the global back out, rather than fighting the export shape.
import "3dmol";
import { fetchAlphaFoldPrediction, plddtColor } from "../utils/alphafold.js";

export default function ProteinStructureViewer({ uniprotId, targetName }) {
  const containerRef = useRef(null);
  const viewerRef = useRef(null);
  const [prediction, setPrediction] = useState(null);
  const [pdbData, setPdbData] = useState(null);
  const [status, setStatus] = useState("loading"); // loading | found | not_found | error

  // Fetch metadata + structure file for this target.
  useEffect(() => {
    let cancelled = false;
    setStatus("loading");
    setPrediction(null);
    setPdbData(null);

    fetchAlphaFoldPrediction(uniprotId)
      .then(async (pred) => {
        if (cancelled) return;
        if (!pred) {
          setStatus("not_found");
          return;
        }
        const pdbRes = await fetch(pred.pdbUrl);
        const data = await pdbRes.text();
        if (cancelled) return;
        setPrediction(pred);
        setPdbData(data);
        setStatus("found");
      })
      .catch(() => {
        if (!cancelled) setStatus("error");
      });

    return () => {
      cancelled = true;
    };
  }, [uniprotId]);

  // Render into the container *after* React has committed it to the DOM —
  // doing this synchronously in the fetch callback above was a bug: state
  // updates are async, so containerRef.current was still null (the
  // "loading" JSX, which has no container div, was still on screen) at the
  // moment the render call actually ran, and it silently no-opped.
  useEffect(() => {
    if (status !== "found" || !pdbData || !containerRef.current) return;

    containerRef.current.innerHTML = "";
    const viewer = window.$3Dmol.createViewer(containerRef.current, {
      backgroundColor: "#111418",
    });
    viewerRef.current = viewer;
    viewer.addModel(pdbData, "pdb");
    viewer.setStyle({}, { cartoon: { colorfunc: (atom) => plddtColor(atom.b) } });
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
        Checking AlphaFold DB for {targetName}…
      </div>
    );
  }

  if (status === "not_found") {
    return (
      <div className="text-sm text-textsecondary py-8 text-center">
        No AlphaFold prediction found for {targetName} ({uniprotId}).
      </div>
    );
  }

  if (status === "error") {
    return (
      <div className="text-sm text-danger py-8 text-center">
        Couldn't reach AlphaFold DB — try again later.
      </div>
    );
  }

  return (
    <div>
      <div
        ref={containerRef}
        style={{ width: "100%", height: 320, position: "relative" }}
      />
      <div className="flex items-center justify-between mt-2 text-xs text-textsecondary">
        <span>
          {prediction.uniprotDescription} — {prediction.organismScientificName}
          {" · "}
          <span className="font-mono">{prediction.sequenceEnd} aa</span>
          {" · "}
          global confidence{" "}
          <span className="font-mono">{prediction.globalMetricValue.toFixed(1)}</span>
        </span>
        <a
          href={`https://alphafold.ebi.ac.uk/entry/${uniprotId}`}
          target="_blank"
          rel="noreferrer"
          className="text-accent hover:underline shrink-0 ml-3"
        >
          View on AlphaFold DB ↗
        </a>
      </div>
      <div className="flex items-center gap-3 mt-2 text-[11px] text-textsecondary">
        <span>Confidence (pLDDT):</span>
        <Swatch color="#0053D6" label="Very high (>90)" />
        <Swatch color="#65CBF3" label="Confident (>70)" />
        <Swatch color="#FFDB13" label="Low (>50)" />
        <Swatch color="#FF7D45" label="Very low" />
      </div>
    </div>
  );
}

function Swatch({ color, label }) {
  return (
    <span className="flex items-center gap-1">
      <span
        className="inline-block w-2.5 h-2.5"
        style={{ backgroundColor: color }}
      />
      {label}
    </span>
  );
}
