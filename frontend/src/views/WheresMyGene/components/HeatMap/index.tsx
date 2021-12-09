import { interpolateYlOrRd } from "d3-scale-chromatic";
import * as echarts from "echarts";
import debounce from "lodash/debounce";
import { useEffect, useMemo, useRef, useState } from "react";
import { EMPTY_ARRAY } from "src/common/constants/utils";
import { CellTypeAndGenes, Gene } from "../../common/types";
import { Container } from "./style";

const DEBOUNCE_MS = 500;

interface Props {
  cellTypes: CellTypeAndGenes[];
  data: CellTypeAndGenes[];
  genes: Gene[];
}

const ELEMENT_ID = "heat-map";

let isChartInitialized = false;
let isObserverInitialized = false;

const MAX_FIRST_PART_LENGTH_PX = 16;

const HEAT_MAP_BASE_WIDTH_PX = 500;
const HEAT_MAP_BASE_HEIGHT_PX = 300;
const HEAT_MAP_BASE_CELL_PX = 20;

export default function HeatMap({
  cellTypes,
  data,
  genes,
}: Props): JSX.Element {
  const [chart, setChart] = useState<echarts.ECharts | null>(null);
  const [xAxisChart, setXAxisChart] = useState<echarts.ECharts | null>(null);
  const [isEchartGLAvailable, setIsEchartGLAvailable] = useState(false);
  const [heatmapWidth, setHeatmapWidth] = useState(getHeatmapWidth(genes));
  const [heatmapHeight, setHeatmapHeight] = useState(
    getHeatmapHeight(cellTypes)
  );

  const [chartData, setChartData] = useState<ChartFormat[]>(EMPTY_ARRAY);

  const [isBottom, setIsBottom] = useState(false);

  const debouncedDataToChartFormat = useMemo(() => {
    return debounce((data, cellTypes, genes) => {
      setChartData(dataToChartFormat(data, cellTypes, genes));
    }, DEBOUNCE_MS);
  }, []);

  useEffect(() => {
    debouncedDataToChartFormat(data, cellTypes, genes);
  }, [data, cellTypes, genes, debouncedDataToChartFormat]);

  // Calculate the min and max of the mean expression of each gene -- end

  const ref = useRef(null);
  const xAxisRef = useRef(null);
  const sentinelRef = useRef(null);

  // Update heatmap size
  useEffect(() => {
    setHeatmapWidth(getHeatmapWidth(genes));
    setHeatmapHeight(getHeatmapHeight(cellTypes));
  }, [cellTypes, genes]);

  // Import echarts-gl
  useEffect(() => {
    importEchartsGl();

    async function importEchartsGl(): Promise<void> {
      await import("echarts-gl");
      setIsEchartGLAvailable(true);
    }
  }, []);

  // set up intersection observer
  useEffect(() => {
    const { current } = sentinelRef;

    if (!current || isObserverInitialized) return;

    const observer = new IntersectionObserver((entries) => {
      if (entries.some((entry) => entry.isIntersecting)) {
        setIsBottom(true);
      } else {
        setIsBottom(false);
      }
    });

    observer.observe(current);

    isObserverInitialized = true;

    return () => observer.disconnect();
  }, [sentinelRef]);

  useEffect(() => {
    const { current } = ref;
    const { current: xAxisCurrent } = xAxisRef;

    if (!current || !xAxisCurrent || isChartInitialized || !isEchartGLAvailable)
      return;

    isChartInitialized = true;
    setChart(echarts.init(current));
    setXAxisChart(echarts.init(xAxisCurrent));
  }, [ref, xAxisRef, isEchartGLAvailable]);

  useEffect(() => {
    if (!chart || !xAxisChart || !isEchartGLAvailable) return;

    const allGeneNames = genes.map((gene) => gene.name);

    const commonSeries = {
      encode: {
        x: "geneIndex",
        y: "cellTypeIndex",
      },
      name: "wmg",
      type: "scatter",
    };

    const commonOptions = {
      dataset: {
        source: chartData,
      },
      large: true,
    };

    chart.setOption({
      ...commonOptions,
      grid: {
        bottom: "100px",
        containLabel: true,
        left: "100px",
        top: "0",
      },
      series: [
        {
          ...commonSeries,
          itemStyle: {
            color(props: { data: { scaledMeanExpression: number } }) {
              const { scaledMeanExpression } = props.data;

              return interpolateYlOrRd(scaledMeanExpression);
            },
          },
          symbolSize: function (props: { percentage: number }) {
            const { percentage } = props;

            return Math.round(MAX_FIRST_PART_LENGTH_PX * percentage);
          },
        },
      ],
      tooltip: {
        formatter(props: { name: string; data: ChartFormat }) {
          const { name, data } = props;
          const {
            geneIndex,
            percentage,
            meanExpression,
            scaledMeanExpression,
          } = data;

          return `
            cell type: ${name}
            <br />
            gene: ${genes[geneIndex].name}
            <br />
            percentage: ${percentage}
            <br />
            mean expression: ${meanExpression}
            <br />
            scaledMeanExpression: ${scaledMeanExpression}
          `;
        },
        position: "top",
      },
      xAxis: [
        {
          axisLabel: { fontSize: isBottom ? 12 : 0, rotate: 90 },
          axisLine: {
            show: false,
          },
          boundaryGap: false,
          data: allGeneNames,
          splitLine: {
            show: true,
          },
          type: "category",
        },
      ],
      yAxis: [
        {
          axisLabel: { rotate: 20 },
          axisLine: {
            show: true,
          },
          data: cellTypes.map((cellType) => cellType.name).reverse(),
          splitLine: {
            show: true,
          },
        },
      ],
    });

    xAxisChart.setOption({
      ...commonOptions,
      grid: {
        containLabel: true,
        left: "100px",
        top: "210px",
      },
      series: [
        {
          ...commonSeries,
          symbolSize: 0,
        },
      ],
      xAxis: [
        {
          axisLabel: { rotate: 90 },
          axisLine: {
            show: false,
          },
          boundaryGap: false,
          data: allGeneNames,
          show: isBottom ? false : true,
          splitLine: {
            show: true,
          },
          type: "category",
        },
      ],
      yAxis: [
        {
          axisLabel: { rotate: 20 },
          axisLine: {
            show: true,
          },
          data: cellTypes.map((cellType) => cellType.name).reverse(),
          show: false,
          splitLine: {
            show: true,
          },
        },
      ],
    });
  }, [
    chart,
    xAxisChart,
    isEchartGLAvailable,
    cellTypes,
    data,
    genes,
    chartData,
    isBottom,
  ]);

  // Resize heatmap
  useEffect(() => {
    chart?.resize();
    xAxisChart?.resize();
  }, [chart, xAxisChart, heatmapHeight, heatmapWidth]);

  return (
    <Container>
      <div
        style={{
          height: "210px",
          position: "sticky",
          top: "0",
          width: heatmapWidth + "px",
          zIndex: 1,
        }}
        ref={xAxisRef}
      />
      <div
        style={{
          height: heatmapHeight + "px",
          width: heatmapWidth + "px",
        }}
        ref={ref}
        id={ELEMENT_ID}
      />
      <div
        style={{
          height: "1px",
          position: "absolute",
          top: `calc(${heatmapHeight}px - 5rem + 180px)`,
          width: heatmapWidth + "px",
        }}
        ref={sentinelRef}
      />
    </Container>
  );
}

interface ChartFormat {
  cellTypeIndex: number;
  geneIndex: number;
  percentage: number;
  meanExpression: number;
  scaledMeanExpression: number;
}

function dataToChartFormat(
  data: CellTypeAndGenes[],
  cellTypes: CellTypeAndGenes[],
  genes: Gene[]
): ChartFormat[] {
  let min = Infinity;
  let max = -Infinity;

  for (const dataPoint of Object.values(data)) {
    for (const expression of Object.values(dataPoint.expressions)) {
      const { me } = expression;

      min = Math.min(min, me);
      max = Math.max(max, me);
    }
  }

  const oldRange = max - min;

  return data.flatMap((dataPoint) => {
    return toChartFormat(dataPoint);
  });

  function toChartFormat(dataPoint: CellTypeAndGenes): ChartFormat[] {
    if (!dataPoint.expressions) return [];

    return Object.entries(dataPoint.expressions).map(
      ([geneName, expression]) => {
        const { pc, me } = expression;

        const scaledMe = (me - min) / oldRange;

        const geneIndex = genes.findIndex((gene) => gene.name === geneName);

        const cellTypeIndex = cellTypes.findIndex(
          (cellType) => cellType.id === dataPoint.id
        );

        return {
          cellTypeIndex,
          geneIndex,
          meanExpression: me,
          percentage: pc,
          scaledMeanExpression: scaledMe,
        } as ChartFormat;
      }
    );
  }
}

function getHeatmapWidth(genes: Gene[]): number {
  return HEAT_MAP_BASE_WIDTH_PX + HEAT_MAP_BASE_CELL_PX * genes.length;
}

function getHeatmapHeight(cellTypes: CellTypeAndGenes[]): number {
  return HEAT_MAP_BASE_HEIGHT_PX + HEAT_MAP_BASE_CELL_PX * cellTypes.length;
}
