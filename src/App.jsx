import { useState } from "react";
import Layout from "./components/Layout.jsx";
import Explorer from "./components/Explorer.jsx";
import Analytics from "./components/Analytics.jsx";
import Analyzer from "./components/Analyzer.jsx";
import TargetLookup from "./components/TargetLookup.jsx";
import NetworkGraph from "./components/NetworkGraph.jsx";
import aptamerData from "./data/aptascope_merged.json";
import datasetStats from "./data/dataset_stats.json";

const TABS = ["Explorer", "Analytics", "Analyzer", "Targets", "Network"];

export default function App() {
  const [activeTab, setActiveTab] = useState("Explorer");

  return (
    <Layout
      tabs={TABS}
      activeTab={activeTab}
      onTabChange={setActiveTab}
      stats={datasetStats}
    >
      {activeTab === "Explorer" && <Explorer data={aptamerData} />}
      {activeTab === "Analytics" && (
        <Analytics data={aptamerData} stats={datasetStats} />
      )}
      {activeTab === "Analyzer" && <Analyzer data={aptamerData} />}
      {activeTab === "Targets" && <TargetLookup data={aptamerData} />}
      {activeTab === "Network" && <NetworkGraph data={aptamerData} />}
    </Layout>
  );
}
