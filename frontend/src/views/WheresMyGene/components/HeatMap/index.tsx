import { interpolateYlOrRd } from "d3-scale-chromatic";
import * as echarts from "echarts";
import { useEffect, useRef, useState } from "react";
import { CellTypeAndGenes, Gene } from "../../common/types";
import { Container } from "./style";

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

export default function HeatMap({
  cellTypes,
  data,
  genes,
}: Props): JSX.Element {
  const [chart, setChart] = useState<echarts.ECharts | null>(null);
  const [isEchartGLAvailable, setIsEchartGLAvailable] = useState(false);
  const [heatmapWidth, setHeatmapWidth] = useState(getHeatmapWidth(genes));
  const [heatmapHeight, setHeatmapHeight] = useState(
    getHeatmapHeight(cellTypes)
  );

  const ref = useRef(null);

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

    if (!current || isChartInitialized || !isEchartGLAvailable) return;

    isChartInitialized = true;
    setChart(echarts.init(current));
  }, [ref, isEchartGLAvailable]);

  useEffect(() => {
    if (!chart || !isEchartGLAvailable) return;

    const allGeneNames = genes.map((gene) => gene.name);

    chart.setOption({
      dataset: {
        source: dataToChartFormat(data, cellTypes, genes),
      },
      grid: {
        bottom: "100px",
        containLabel: true,
        left: "100px",
        top: "0",
      },
      large: true,
      series: [
        {
          encode: {
            x: "geneIndex",
            y: "cellTypeIndex",
          },
          itemStyle: {
            color(props: { data: { meanExpression: number } }) {
              const { meanExpression } = props.data;

              return interpolateYlOrRd(meanExpression);
            },
          },
          name: "wmg",
          symbolSize: function (props: { percentage: number }) {
            const { percentage } = props;

            return Math.round(MAX_FIRST_PART_LENGTH_PX * percentage);
          },
          type: "scatter",
        },
      ],
      tooltip: {
        formatter(props: { name: string; data: ChartFormat }) {
          const { name, data } = props;
          const { geneIndex, percentage, meanExpression } = data;

          return `
            cell type: ${name}
            <br />
            gene: ${genes[geneIndex].name}
            <br />
            percentage: ${percentage}
            <br />
            mean expression: ${meanExpression}
          `;
        },
        position: "top",
      },
      xAxis: [
        {
          axisLabel: { rotate: 90 },
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
  }, [chart, isEchartGLAvailable, cellTypes, data, genes]);

  // Resize heatmap
  useEffect(() => {
    chart?.resize();
  }, [chart, heatmapHeight, heatmapWidth]);

  return (
    <Container>
      <div
        style={{
          height: heatmapHeight + "px",
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
}

function dataToChartFormat(
  data: CellTypeAndGenes[],
  cellTypes: CellTypeAndGenes[],
  genes: Gene[]
): ChartFormat[] {
  return data.flatMap((dataPoint) => {
    return toChartFormat(dataPoint);
  });

  function toChartFormat(dataPoint: CellTypeAndGenes): ChartFormat[] {
    if (!dataPoint.expressions) return [];

    return Object.entries(dataPoint.expressions).map(
      ([geneName, expression]) => {
        const { pc, me } = expression;

        const geneIndex = genes.findIndex((gene) => gene.name === geneName);

        const cellTypeIndex = cellTypes.findIndex(
          (cellType) => cellType.id === dataPoint.id
        );

        return {
          cellTypeIndex,
          geneIndex,
          meanExpression: me,
          percentage: pc,
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
