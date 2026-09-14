import { useEffect, useState } from "react";
import Layout from "./components/Layout.jsx";
import Explorer from "./components/Explorer.jsx";
import Analytics from "./components/Analytics.jsx";
import Analyzer from "./components/Analyzer.jsx";
import TargetLookup from "./components/TargetLookup.jsx";
import NetworkGraph from "./components/NetworkGraph.jsx";
import datasetStats from "./data/dataset_stats.json";

const TABS = ["Explorer", "Analytics", "Analyzer", "Targets", "Network"];

// Served from public/ as a plain static asset (16MB+) instead of bundled
// into the JS — keeps the parsed/executed JS payload small and lets the
// browser cache this file independently of app-code deploys.
const DATA_URL = `${import.meta.env.BASE_URL}data/aptascope_merged.json`;

export default function App() {
  const [activeTab, setActiveTab] = useState("Explorer");
  const [aptamerData, setAptamerData] = useState(null);
  const [loadError, setLoadError] = useState(false);

  useEffect(() => {
    fetch(DATA_URL)
      .then((res) => {
        if (!res.ok) throw new Error(`${res.status} ${res.statusText}`);
        return res.json();
      })
      .then(setAptamerData)
      .catch(() => setLoadError(true));
  }, []);

  return (
    <Layout tabs={TABS} activeTab={activeTab} onTabChange={setActiveTab} stats={datasetStats}>
      {!aptamerData ? (
        <div className="py-24 text-center text-sm text-textsecondary">
          {loadError ? "Failed to load the aptamer dataset. Try reloading." : "Loading dataset…"}
        </div>
      ) : (
        <>
          {activeTab === "Explorer" && <Explorer data={aptamerData} />}
          {activeTab === "Analytics" && (
            <Analytics data={aptamerData} stats={datasetStats} />
          )}
          {activeTab === "Analyzer" && <Analyzer data={aptamerData} />}
          {activeTab === "Targets" && <TargetLookup data={aptamerData} />}
          {activeTab === "Network" && <NetworkGraph data={aptamerData} />}
        </>
      )}
    </Layout>
  );
}
