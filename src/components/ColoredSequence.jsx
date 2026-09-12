import { useState } from "react";
import { nucleotideColor } from "../utils/sequence.js";

export function CopyButton({ text, className = "" }) {
  const [copied, setCopied] = useState(false);

  const handleCopy = async (e) => {
    e.stopPropagation();
    try {
      await navigator.clipboard.writeText(text);
      setCopied(true);
      setTimeout(() => setCopied(false), 1500);
    } catch {
      // clipboard API unavailable — silently ignore
    }
  };

  return (
    <button
      onClick={handleCopy}
      className={`text-xs px-2 py-0.5 rounded border border-border text-textsecondary hover:text-textprimary hover:border-accent transition-colors ${className}`}
      title="Copy sequence"
    >
      {copied ? "Copied" : "Copy"}
    </button>
  );
}

export default function ColoredSequence({ sequence, className = "" }) {
  return (
    <span className={`font-mono tracking-[0.05em] break-all ${className}`}>
      {sequence.split("").map((base, i) => (
        <span key={i} style={{ color: nucleotideColor(base) }}>
          {base}
        </span>
      ))}
    </span>
  );
}
