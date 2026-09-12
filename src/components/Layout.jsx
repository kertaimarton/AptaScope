import { useState } from "react";
import Logo from "./Logo.jsx";

export default function Layout({ tabs, activeTab, onTabChange, stats, children }) {
  const [showAbout, setShowAbout] = useState(false);

  return (
    <div className="min-h-screen bg-bg text-textprimary">
      <header className="border-b border-border">
        <div className="max-w-[1280px] mx-auto px-6">
          <div className="flex items-center justify-between h-16">
            <div className="flex items-center gap-3">
              <Logo />
              <span className="hidden sm:inline text-textsecondary text-xs uppercase tracking-wider border border-border px-3 py-0.5 font-mono">
                {stats.total_records.toLocaleString()} apt · {stats.unique_targets.toLocaleString()} tgt
              </span>
            </div>
            <nav className="flex items-center gap-6">
              {tabs.map((tab) => (
                <button
                  key={tab}
                  onClick={() => onTabChange(tab)}
                  className={`relative py-1.5 text-sm uppercase tracking-wider transition-colors ${
                    activeTab === tab
                      ? "text-textprimary"
                      : "text-textsecondary hover:text-textprimary"
                  }`}
                >
                  {tab}
                  {activeTab === tab && (
                    <span className="absolute left-0 right-0 -bottom-[1px] h-[2px] bg-accent" />
                  )}
                </button>
              ))}
            </nav>
          </div>
        </div>
      </header>

      <main className="max-w-[1280px] mx-auto px-6 py-8">{children}</main>

      <footer className="border-t border-border mt-16">
        <div className="max-w-[1280px] mx-auto px-6 py-6 flex items-center justify-between text-xs text-textsecondary">
          <span>
            AptaScope — research exploration tool, not a clinical or regulatory instrument.
          </span>
          <button
            onClick={() => setShowAbout(true)}
            className="text-xs uppercase tracking-wider border border-border px-3 py-1 text-textsecondary hover:text-textprimary hover:border-accent transition-colors"
          >
            About &amp; Data Sources
          </button>
        </div>
      </footer>

      {showAbout && <AboutModal onClose={() => setShowAbout(false)} />}
    </div>
  );
}

function AboutModal({ onClose }) {
  return (
    <div
      className="fixed inset-0 bg-black/60 flex items-center justify-center z-50 p-4"
      onClick={onClose}
    >
      <div
        className="bg-surface border border-border max-w-xl w-full p-6 space-y-4"
        onClick={(e) => e.stopPropagation()}
      >
        <div className="flex items-center justify-between">
          <div>
            <h2 className="text-lg font-semibold">About AptaScope</h2>
            <div className="h-[2px] w-8 bg-accent mt-1.5" />
          </div>
          <button
            onClick={onClose}
            className="text-textsecondary hover:text-textprimary text-xl leading-none"
          >
            ×
          </button>
        </div>

        <p className="text-[13px] text-textsecondary leading-relaxed">
          AptaScope is an in-browser tool for analyzing aptamer sequences —
          short DNA/RNA molecules that bind specific targets — against a
          database of experimentally validated aptamer-target interactions.
          It computes biophysical features for any pasted sequence, compares
          it to known high-affinity binders, and surfaces the most similar
          characterized aptamers.
        </p>

        <div className="text-sm space-y-1.5">
          <p className="text-xs uppercase tracking-wider text-textsecondary">Data sources</p>
          <ul className="text-textsecondary space-y-1 list-disc list-inside">
            <li>
              <a
                href="https://zenodo.org/records/8264921"
                target="_blank"
                rel="noreferrer"
                className="text-accent hover:underline"
              >
                UTexas Aptamer Database
              </a>{" "}
              (Zenodo)
            </li>
            <li>
              <a
                href="https://aptamer.ribocentre.org/"
                target="_blank"
                rel="noreferrer"
                className="text-accent hover:underline"
              >
                Ribocentre Aptamer
              </a>
            </li>
          </ul>
        </div>

        <p className="text-[13px] text-textsecondary">
          <span className="font-medium text-textprimary">Kd (dissociation constant)</span>{" "}
          measures binding affinity — lower values mean tighter, stronger binding.
        </p>

        <p className="text-xs uppercase tracking-wider text-textsecondary border-t border-border pt-3">
          This tool is for research exploration · validate computationally-derived
          insights experimentally
        </p>
      </div>
    </div>
  );
}
