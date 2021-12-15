import { interpolateYlOrRd } from "d3-scale-chromatic";
import * as echarts from "echarts";
import debounce from "lodash/debounce";
import { useEffect, useMemo, useRef, useState } from "react";
import { EMPTY_OBJECT } from "src/common/constants/utils";
import { CellTypeAndGenes, Gene } from "../../common/types";
import { Container } from "./style";

interface ChartProps {
  chartData: ChartFormat[];
  geneNames: string[];
  cellTypeNames: string[];
}

const DEBOUNCE_MS = 2 * 1000;

interface Props {
  cellTypes: CellTypeAndGenes[];
  data: CellTypeAndGenes[];
  genes: Gene[];
}

const ELEMENT_ID = "heat-map";

let isChartInitialized = false;

const MAX_FIRST_PART_LENGTH_PX = 16;

const HEAT_MAP_BASE_WIDTH_PX = 500;
const HEAT_MAP_BASE_HEIGHT_PX = 300;
const HEAT_MAP_BASE_CELL_PX = 20;

const COMMON_SERIES = {
  emphasis: { itemStyle: { color: "inherit" }, scale: false },
  encode: {
    x: "geneIndex",
    y: "cellTypeIndex",
  },
  legendHoverLink: false,
  name: "wmg",
  type: "scatter",
};

const COMMON_AXIS_POINTER = {
  link: [
    {
      xAxisIndex: "all",
    },
  ],
  show: true,
  triggerOn: "click",
  triggerTooltip: false,
};

export default function HeatMap({
  cellTypes,
  data,
  genes,
}: Props): JSX.Element {
  const [chart, setChart] = useState<echarts.ECharts | null>(null);
  const [xAxisChart, setXAxisChart] = useState<echarts.ECharts | null>(null);
  const [yAxisChart, setYAxisChart] = useState<echarts.ECharts | null>(null);
  const [isEchartGLAvailable, setIsEchartGLAvailable] = useState(false);
  const [heatmapWidth, setHeatmapWidth] = useState(getHeatmapWidth(genes));
  const [heatmapHeight, setHeatmapHeight] = useState(
    getHeatmapHeight(cellTypes)
  );

  const [chartProps, setChartProps] = useState<ChartProps | null>(null);

  // Loading state
  const [isLoading, setIsLoading] = useState(false);

  useEffect(() => {
    setIsLoading(true);
  }, [genes]);

  useEffect(() => {
    if (!chart) return;

    if (isLoading) {
      chart.showLoading();

      // DEBUG
      // DEBUG
      // DEBUG
      console.log("--show loading");
    } else {
      // DEBUG
      // DEBUG
      // DEBUG
      console.log("--hide loading");

      chart.hideLoading();
    }
  }, [chart, isLoading]);

  const debouncedDataToChartFormat = useMemo(() => {
    return debounce(
      (data, cellTypes, genes) => {
        // DEBUG
        // DEBUG
        // DEBUG
        console.log("------------SETTING!!");

        setChartProps({
          cellTypeNames: getCellTypeNames(cellTypes),
          chartData: dataToChartFormat(data, cellTypes, genes),
          geneNames: getGeneNames(genes),
        });

        setIsLoading(false);
      },
      DEBOUNCE_MS,
      { leading: false }
    );
  }, []);

  useEffect(() => {
    debouncedDataToChartFormat(data, cellTypes, genes);
  }, [data, cellTypes, genes, debouncedDataToChartFormat]);

  // Calculate the min and max of the mean expression of each gene -- end

  const ref = useRef(null);
  const xAxisRef = useRef(null);
  const yAxisRef = useRef(null);

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

  useEffect(() => {
    const { current } = ref;
    const { current: xAxisCurrent } = xAxisRef;
    const { current: yAxisCurrent } = yAxisRef;

    if (
      !current ||
      !xAxisCurrent ||
      !yAxisCurrent ||
      isChartInitialized ||
      !isEchartGLAvailable
    ) {
      return;
    }

    isChartInitialized = true;

    const chart = echarts.init(current, EMPTY_OBJECT, { useDirtyRect: true });
    const xAxisChart = echarts.init(xAxisCurrent, EMPTY_OBJECT, {
      useDirtyRect: true,
    });
    const yAxisChart = echarts.init(yAxisCurrent, EMPTY_OBJECT, {
      useDirtyRect: true,
    });

    echarts.connect([chart, xAxisChart, yAxisChart]);

    setChart(chart);
    setXAxisChart(xAxisChart);
    setYAxisChart(yAxisChart);
  }, [ref, xAxisRef, isEchartGLAvailable]);

  useEffect(() => {
    if (
      !chart ||
      !chartProps ||
      !xAxisChart ||
      !yAxisChart ||
      !isEchartGLAvailable
    ) {
      return;
    }

    const { chartData, cellTypeNames, geneNames } = chartProps;

    const commonOptions = {
      animation: false,
      dataset: {
        source: chartData,
      },
      hoverLayerThreshold: 10,
      progressiveThreshold: 2000,
    };

    // DEBUG
    // DEBUG
    // DEBUG
    console.log("----- rendering...");
    console.log("----- chartData", chartData);
    console.log("----- geneNames", geneNames);
    console.log("----- cellTypeNames", cellTypeNames);

    chart.setOption({
      ...commonOptions,
      axisPointer: { ...COMMON_AXIS_POINTER, label: { show: false } },
      grid: {
        bottom: "100px",
        left: "300px",
        top: "200px",
      },
      series: [
        {
          ...COMMON_SERIES,
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

          if (!data) return false;

          const {
            geneIndex,
            percentage,
            meanExpression,
            scaledMeanExpression,
          } = data;

          return `
            cell type: ${name}
            <br />
            gene: ${geneNames[geneIndex]}
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
          axisLabel: { fontSize: 0, rotate: 90 },
          axisLine: {
            show: false,
          },
          axisTick: {
            show: false,
          },
          boundaryGap: false,
          data: geneNames,
          splitLine: {
            show: true,
          },
          type: "category",
        },
      ],
      yAxis: [
        {
          axisLabel: { fontSize: 0, rotate: 20 },
          axisLine: {
            show: false,
          },
          axisTick: {
            show: false,
          },
          boundaryGap: false,
          data: cellTypeNames,
          splitLine: {
            show: true,
          },
        },
      ],
    });

    xAxisChart.setOption({
      ...commonOptions,
      axisPointer: { ...COMMON_AXIS_POINTER },
      grid: {
        bottom: "0",
        left: "300px",
        top: "200px",
      },
      series: [
        {
          ...COMMON_SERIES,
          symbolSize: 0,
        },
      ],
      xAxis: [
        {
          axisLabel: {
            overflow: "truncate",
            rotate: 270,
            verticalAlign: "bottom",
            width: 200,
          },
          boundaryGap: false,
          data: geneNames,
          position: "top",
          type: "category",
        },
      ],
      yAxis: [
        {
          axisLabel: { fontSize: 0, rotate: 20 },
          axisLine: {
            show: false,
          },
          axisTick: {
            show: false,
          },
          boundaryGap: false,
          data: cellTypeNames,
          splitLine: {
            show: false,
          },
        },
      ],
    });

    yAxisChart.setOption({
      ...commonOptions,
      axisPointer: { ...COMMON_AXIS_POINTER },
      grid: {
        bottom: "300px",
        left: "300px",
        right: 0,
        top: 0,
      },
      series: [
        {
          ...COMMON_SERIES,
          symbolSize: 0,
        },
      ],
      xAxis: [
        {
          axisLabel: { fontSize: 0, rotate: 90 },
          axisLine: {
            show: false,
          },
          axisTick: {
            show: false,
          },
          boundaryGap: false,
          data: geneNames,
          splitLine: {
            show: false,
          },
          type: "category",
        },
      ],
      yAxis: [
        {
          axisLabel: { align: "right", rotate: 20 },
          axisLine: {
            show: false,
          },
          axisTick: {
            alignWithLabel: true,
          },
          boundaryGap: false,
          data: cellTypeNames,
          splitLine: {
            show: false,
          },
        },
      ],
    });

    chart.resize();
    xAxisChart.resize();
    yAxisChart.resize();
  }, [chart, xAxisChart, yAxisChart, isEchartGLAvailable, chartProps]);

  return (
    <Container>
      <div
        style={{
          backgroundColor: "rgba(255, 255, 255, 0.8)",
          height: "200px",
          position: "sticky",
          top: "0",
          width: heatmapWidth + "px",
          zIndex: 2,
        }}
        ref={xAxisRef}
      />
      <div
        style={{
          backgroundColor: "rgba(255, 255, 255, 0.8)",
          height: heatmapHeight + "px",
          left: 0,
          position: "sticky",
          top: 0,
          width: "300px",
          zIndex: 1,
        }}
        ref={yAxisRef}
      />
      <div
        style={{
          height: heatmapHeight + "px",
          left: "0",
          position: "absolute",
          top: "0",
          width: heatmapWidth + "px",
        }}
        ref={ref}
        id={ELEMENT_ID}
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

function getCellTypeNames(cellTypes: CellTypeAndGenes[]): string[] {
  return cellTypes.map((cellType) => cellType.name).reverse();
}

function getGeneNames(genes: Gene[]): string[] {
  return genes.map((gene) => gene.name);
}
