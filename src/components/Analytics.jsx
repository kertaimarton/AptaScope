import { useMemo } from "react";
import {
  ResponsiveContainer,
  ScatterChart,
  Scatter,
  XAxis,
  YAxis,
  ZAxis,
  CartesianGrid,
  Tooltip,
  Legend,
  Line,
  ComposedChart,
  BarChart,
  Bar,
  Cell,
  ReferenceLine,
  LineChart,
} from "recharts";
import { median, quantile, movingAverage } from "../utils/stats.js";
import PageHeading from "./PageHeading.jsx";

const COLORS = {
  dna: "#58A6FF",
  rna: "#F0883E",
  accent: "#E0576B",
  success: "#238551",
  warning: "#C87619",
  danger: "#CD4246",
  grid: "#2F343C",
  text: "#8F99A8",
  nucleotideC: "#58A6FF", // fixed bioinformatics convention (C=blue), independent of the brand accent
};

const AXIS_STYLE = { fontSize: 11, fill: COLORS.text, fontFamily: "JetBrains Mono, monospace" };

function ChartCard({ title, interpretation, children, height = 320 }) {
  return (
    <div className="bg-surface border border-border rounded-md p-4">
      <h3 className="text-xs uppercase tracking-wider font-semibold mb-3">{title}</h3>
      <div style={{ width: "100%", height }}>{children}</div>
      <p className="text-[13px] text-textsecondary mt-3 border-t border-border pt-2">
        {interpretation}
      </p>
    </div>
  );
}

function TooltipBox({ active, payload }) {
  if (!active || !payload || !payload.length) return null;
  return (
    <div className="bg-bg border border-border px-2 py-1.5 text-xs font-mono space-y-0.5">
      {payload.map((p, i) => (
        <div key={i} style={{ color: p.color }}>
          {p.name}: {typeof p.value === "number" ? p.value.toFixed(2) : p.value}
        </div>
      ))}
    </div>
  );
}

export default function Analytics({ data }) {
  const withKd = useMemo(() => data.filter((r) => r.kd_nM !== null && r.kd_nM > 0), [data]);

  // Chart 1: GC content vs binding affinity
  const gcScatter = useMemo(() => {
    const dna = [];
    const rna = [];
    for (const r of withKd) {
      const point = { x: r.gc_content, y: -Math.log10(r.kd_nM), target: r.target_name };
      (r.type === "DNA" ? dna : rna).push(point);
    }
    const trend = movingAverage([...dna, ...rna], 25);
    return { dna, rna, trend };
  }, [withKd]);

  // Chart 2: length vs binding affinity
  const lengthScatter = useMemo(() => {
    return withKd.map((r) => ({ x: r.length, y: -Math.log10(r.kd_nM) }));
  }, [withKd]);

  // Chart 3: affinity distribution histogram
  const affinityHistogram = useMemo(() => {
    const logs = withKd.map((r) => Math.log10(r.kd_nM));
    const min = Math.min(...logs);
    const max = Math.max(...logs);
    const bins = 24;
    const width = (max - min) / bins;
    const counts = new Array(bins).fill(0);
    for (const v of logs) {
      const idx = Math.min(bins - 1, Math.floor((v - min) / width));
      counts[idx] += 1;
    }
    const med = median(logs);
    const sorted = [...logs].sort((a, b) => a - b);
    const q1 = quantile(sorted, 0.25);
    const q3 = quantile(sorted, 0.75);
    return {
      bars: counts.map((count, i) => ({
        log10Kd: min + width * (i + 0.5),
        count,
      })),
      median: med,
      q1,
      q3,
    };
  }, [withKd]);

  // Chart 4: top 15 targets by count, colored by median Kd
  const topTargets = useMemo(() => {
    const byTarget = new Map();
    for (const r of data) {
      if (!byTarget.has(r.target_name)) byTarget.set(r.target_name, []);
      byTarget.get(r.target_name).push(r);
    }
    const rows = [...byTarget.entries()]
      .map(([name, recs]) => {
        const kds = recs.map((r) => r.kd_nM).filter((k) => k !== null);
        return {
          name: name.length > 28 ? name.slice(0, 26) + "…" : name,
          count: recs.length,
          medianKd: kds.length ? median(kds) : null,
        };
      })
      .sort((a, b) => b.count - a.count)
      .slice(0, 15)
      .reverse(); // so the largest ends up at top of horizontal bar chart
    return rows;
  }, [data]);

  function targetBarColor(medianKd) {
    if (medianKd === null) return COLORS.text;
    if (medianKd < 10) return COLORS.success;
    if (medianKd <= 1000) return COLORS.warning;
    return COLORS.danger;
  }

  // Chart 5: nucleotide composition by target type (protein / small_molecule / cell)
  const compositionByType = useMemo(() => {
    const groups = ["protein", "small_molecule", "cell"];
    const labels = { protein: "Protein", small_molecule: "Small Molecule", cell: "Cell" };
    return groups.map((g) => {
      const recs = data.filter((r) => r.target_type === g);
      const avg = (key) =>
        recs.length ? recs.reduce((s, r) => s + r[key], 0) / recs.length : 0;
      return {
        name: labels[g],
        A: avg("a_freq"),
        T: avg("t_freq"),
        G: avg("g_freq"),
        C: avg("c_freq"),
      };
    });
  }, [data]);

  // Chart 6: timeline of aptamers published per year
  const timeline = useMemo(() => {
    const byYear = new Map();
    for (const r of data) {
      if (!r.year) continue;
      byYear.set(r.year, (byYear.get(r.year) || 0) + 1);
    }
    return [...byYear.entries()]
      .sort((a, b) => a[0] - b[0])
      .map(([year, count]) => ({ year, count }));
  }, [data]);

  return (
    <div className="space-y-4">
      <PageHeading>Analytics</PageHeading>

      <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
        <ChartCard
          title="GC Content vs. Binding Affinity"
          interpretation="Points cluster toward the upper-middle when high-affinity binders favor moderate GC content; a flat trend line means GC content alone doesn't predict affinity in this dataset."
        >
          <ResponsiveContainer>
            <ComposedChart margin={{ top: 10, right: 10, bottom: 10, left: 0 }}>
              <CartesianGrid stroke={COLORS.grid} strokeDasharray="3 3" />
              <XAxis
                type="number"
                dataKey="x"
                domain={[0, 1]}
                name="GC Content"
                tick={AXIS_STYLE}
                label={{ value: "GC Content", position: "insideBottom", offset: -5, style: AXIS_STYLE }}
              />
              <YAxis
                type="number"
                dataKey="y"
                name="-log10(Kd)"
                tick={AXIS_STYLE}
                label={{ value: "-log10(Kd nM)", angle: -90, position: "insideLeft", style: AXIS_STYLE }}
              />
              <Tooltip content={<TooltipBox />} cursor={{ strokeDasharray: "3 3" }} />
              <Legend wrapperStyle={{ fontSize: 11 }} />
              <Scatter name="DNA" data={gcScatter.dna} fill={COLORS.dna} opacity={0.6} />
              <Scatter name="RNA" data={gcScatter.rna} fill={COLORS.rna} opacity={0.6} />
              <Line
                type="monotone"
                dataKey="y"
                data={gcScatter.trend}
                stroke={COLORS.text}
                dot={false}
                name="Trend"
                strokeWidth={2}
              />
            </ComposedChart>
          </ResponsiveContainer>
        </ChartCard>

        <ChartCard
          title="Length vs. Binding Affinity"
          interpretation="A visible peak in high-affinity points at a particular length range suggests a structural 'sweet spot' for this dataset's aptamer families."
        >
          <ResponsiveContainer>
            <ScatterChart margin={{ top: 10, right: 10, bottom: 10, left: 0 }}>
              <CartesianGrid stroke={COLORS.grid} strokeDasharray="3 3" />
              <XAxis
                type="number"
                dataKey="x"
                name="Length"
                tick={AXIS_STYLE}
                label={{ value: "Length (nt)", position: "insideBottom", offset: -5, style: AXIS_STYLE }}
              />
              <YAxis
                type="number"
                dataKey="y"
                name="-log10(Kd)"
                tick={AXIS_STYLE}
                label={{ value: "-log10(Kd nM)", angle: -90, position: "insideLeft", style: AXIS_STYLE }}
              />
              <ZAxis range={[20, 20]} />
              <Tooltip content={<TooltipBox />} cursor={{ strokeDasharray: "3 3" }} />
              <Scatter data={lengthScatter} fill={COLORS.accent} opacity={0.5} />
            </ScatterChart>
          </ResponsiveContainer>
        </ChartCard>

        <ChartCard
          title="Affinity Distribution"
          interpretation={`Median log10(Kd) is ${affinityHistogram.median?.toFixed(2)}; the interquartile range (dashed lines) shows where the middle 50% of binders fall.`}
        >
          <ResponsiveContainer>
            <BarChart data={affinityHistogram.bars} margin={{ top: 10, right: 10, bottom: 10, left: 0 }}>
              <CartesianGrid stroke={COLORS.grid} strokeDasharray="3 3" />
              <XAxis
                dataKey="log10Kd"
                tickFormatter={(v) => v.toFixed(1)}
                tick={AXIS_STYLE}
                label={{ value: "log10(Kd nM)", position: "insideBottom", offset: -5, style: AXIS_STYLE }}
              />
              <YAxis tick={AXIS_STYLE} />
              <Tooltip content={<TooltipBox />} />
              <Bar dataKey="count" fill={COLORS.accent} />
              {affinityHistogram.median !== null && (
                <ReferenceLine x={affinityHistogram.median} stroke={COLORS.success} strokeWidth={2} />
              )}
              {affinityHistogram.q1 !== null && (
                <ReferenceLine x={affinityHistogram.q1} stroke={COLORS.text} strokeDasharray="4 4" />
              )}
              {affinityHistogram.q3 !== null && (
                <ReferenceLine x={affinityHistogram.q3} stroke={COLORS.text} strokeDasharray="4 4" />
              )}
            </BarChart>
          </ResponsiveContainer>
        </ChartCard>

        <ChartCard
          title="Top 15 Most-Studied Targets"
          interpretation="Bar color reflects median Kd for that target (green < 10 nM, amber 10–1000 nM, red > 1000 nM) — heavily studied targets aren't always the tightest binders."
        >
          <ResponsiveContainer>
            <BarChart
              data={topTargets}
              layout="vertical"
              margin={{ top: 10, right: 20, bottom: 10, left: 10 }}
            >
              <CartesianGrid stroke={COLORS.grid} strokeDasharray="3 3" horizontal={false} />
              <XAxis type="number" tick={AXIS_STYLE} />
              <YAxis
                type="category"
                dataKey="name"
                width={160}
                tick={{ ...AXIS_STYLE, fontSize: 10 }}
              />
              <Tooltip content={<TooltipBox />} />
              <Bar dataKey="count" name="Aptamers">
                {topTargets.map((t, i) => (
                  <Cell key={i} fill={targetBarColor(t.medianKd)} />
                ))}
              </Bar>
            </BarChart>
          </ResponsiveContainer>
        </ChartCard>

        <ChartCard
          title="Nucleotide Composition by Target Type"
          interpretation="Systematic differences in base composition across target classes hint that different target types favor different aptamer chemistries or motifs."
        >
          <ResponsiveContainer>
            <BarChart data={compositionByType} margin={{ top: 10, right: 10, bottom: 10, left: 0 }}>
              <CartesianGrid stroke={COLORS.grid} strokeDasharray="3 3" />
              <XAxis dataKey="name" tick={AXIS_STYLE} />
              <YAxis tick={AXIS_STYLE} tickFormatter={(v) => `${(v * 100).toFixed(0)}%`} />
              <Tooltip content={<TooltipBox />} />
              <Legend wrapperStyle={{ fontSize: 11 }} />
              <Bar dataKey="A" fill={COLORS.success} />
              <Bar dataKey="T" fill={COLORS.danger} />
              <Bar dataKey="G" fill={COLORS.warning} />
              <Bar dataKey="C" fill={COLORS.nucleotideC} />
            </BarChart>
          </ResponsiveContainer>
        </ChartCard>

        <ChartCard
          title="Aptamers Published per Year"
          interpretation="Publication volume over time reflects both SELEX adoption trends and the coverage of the underlying literature databases, not necessarily field-wide activity."
        >
          <ResponsiveContainer>
            <LineChart data={timeline} margin={{ top: 10, right: 10, bottom: 10, left: 0 }}>
              <CartesianGrid stroke={COLORS.grid} strokeDasharray="3 3" />
              <XAxis dataKey="year" tick={AXIS_STYLE} />
              <YAxis tick={AXIS_STYLE} />
              <Tooltip content={<TooltipBox />} />
              <Line type="monotone" dataKey="count" stroke={COLORS.accent} strokeWidth={2} dot={false} />
            </LineChart>
          </ResponsiveContainer>
        </ChartCard>
      </div>
    </div>
  );
}
