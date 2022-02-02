import { interpolateYlOrRd } from "d3-scale-chromatic";
import {
  DefaultLabelFormatterCallbackParams,
  EChartsOption,
  ScatterSeriesOption,
} from "echarts";
import { State } from "../../common/store";
import { CellTypeSummary } from "../../common/types";

export const MAX_FIRST_PART_LENGTH_PX = 16;
export const X_AXIS_CHART_HEIGHT = "200px";
export const Y_AXIS_CHART_WIDTH = "300px";

interface CreateOptionsProps {
  geneNames: string[];
  cellTypeNames: string[];
  commonOptions: EChartsOption;
}

const COMMON_SERIES: ScatterSeriesOption = {
  emphasis: { itemStyle: { color: "inherit" }, scale: false },
  encode: {
    x: "geneIndex",
    y: "cellTypeIndex",
  },
  legendHoverLink: false,
  name: "wmg",
  type: "scatter",
};

const COMMON_AXIS_POINTER: EChartsOption["axisPointer"] = {
  link: [
    {
      xAxisIndex: "all",
    },
  ],
  show: true,
  triggerOn: "click",
  triggerTooltip: false,
};

export function createChartOptions({
  geneNames,
  cellTypeNames,
  commonOptions,
}: CreateOptionsProps): EChartsOption {
  return {
    ...commonOptions,
    axisPointer: { ...COMMON_AXIS_POINTER, label: { show: false } },
    grid: {
      bottom: "100px",
      left: Y_AXIS_CHART_WIDTH,
      top: X_AXIS_CHART_HEIGHT,
    },
    series: [
      {
        ...COMMON_SERIES,
        itemStyle: {
          color(props: DefaultLabelFormatterCallbackParams) {
            const { scaledMeanExpression } = props.data as {
              scaledMeanExpression: number;
            };

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
      formatter(props) {
        const { name: cellTypeName, data } = props as unknown as {
          name: string;
          data: ChartFormat;
        };

        if (!data) return "";

        const { geneIndex, percentage, meanExpression, scaledMeanExpression } =
          data;

        const { name } = deserializeCellTypeName(cellTypeName);

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
        boundaryGap: true,
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
        boundaryGap: true,
        data: cellTypeNames,
        splitLine: {
          show: true,
        },
      },
    ],
  };
}

interface CreateXAxisOptionsProps extends CreateOptionsProps {
  genesToDelete: string[];
}

export function createXAxisOptions({
  geneNames,
  cellTypeNames,
  commonOptions,
  genesToDelete,
}: CreateXAxisOptionsProps): EChartsOption {
  return {
    ...commonOptions,
    axisPointer: { ...COMMON_AXIS_POINTER },
    grid: {
      bottom: "0",
      left: "300px",
      top: "300px",
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
          formatter(value) {
            return genesToDelete.includes(value)
              ? `{selected|${value}}`
              : value;
          },
          rich: {
            selected: {
              color: "red",
              fontWeight: "bold",
            },
          },
          rotate: 270,
          verticalAlign: "middle",
          width: 200,
        },
        axisTick: {
          alignWithLabel: true,
        },
        boundaryGap: true,
        data: geneNames,
        position: "top",
        triggerEvent: true,
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
        boundaryGap: true,
        data: cellTypeNames,
        splitLine: {
          show: false,
        },
      },
    ],
  };
}

interface CreateYAxisOptionsProps extends CreateOptionsProps {
  cellTypeIdsToDelete: string[];
}

export function createYAxisOptions({
  geneNames,
  cellTypeNames,
  commonOptions,
  cellTypeIdsToDelete,
}: CreateYAxisOptionsProps): EChartsOption {
  return {
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
        boundaryGap: true,
        data: geneNames,
        splitLine: {
          show: false,
        },
        type: "category",
      },
    ],
    yAxis: [
      {
        axisLabel: {
          align: "right",
          formatter(value) {
            const [id, , name] = value.split("~");

            return cellTypeIdsToDelete.includes(id)
              ? `{selected|${name}}`
              : name;
          },
          rich: {
            selected: {
              color: "red",
              fontWeight: "bold",
            },
          },
          rotate: 20,
        },
        axisLine: {
          show: false,
        },
        axisTick: {
          alignWithLabel: true,
        },
        boundaryGap: true,
        data: cellTypeNames,
        splitLine: {
          show: false,
        },
        triggerEvent: true,
      },
    ],
  };
}

export interface ChartFormat {
  cellTypeIndex: number;
  geneIndex: number;
  percentage: number;
  meanExpression: number;
  scaledMeanExpression: number;
}

export function dataToChartFormat(
  data: CellTypeSummary[],
  cellTypes: CellTypeSummary[],
  genes: State["selectedGenes"]
): ChartFormat[] {
  let min = Infinity;
  let max = -Infinity;

  for (const dataPoint of Object.values(data)) {
    const { geneExpressions } = dataPoint;

    if (!geneExpressions) continue;

    for (const geneExpression of Object.values(geneExpressions)) {
      const { meanExpression } = geneExpression;

      min = Math.min(min, meanExpression);
      max = Math.max(max, meanExpression);
    }
  }

  const oldRange = max - min;

  const result = data.flatMap((dataPoint) => {
    return toChartFormat(dataPoint);
  });

  return result;

  function toChartFormat(dataPoint: CellTypeSummary): ChartFormat[] {
    const { geneExpressions } = dataPoint;

    if (!geneExpressions) return [];

    return Object.entries(geneExpressions).map(([geneName, geneExpression]) => {
      const { percentage, meanExpression } = geneExpression;

      const scaledMeanExpression = (meanExpression - min) / oldRange;

      const geneIndex = genes.findIndex((gene) => gene === geneName);

      const cellTypeIndex = cellTypes.findIndex(
        (cellType) => cellType.id === dataPoint.id
      );

      return {
        cellTypeIndex,
        geneIndex,
        meanExpression,
        percentage,
        scaledMeanExpression,
      } as ChartFormat;
    });
  }
}

const HEAT_MAP_BASE_WIDTH_PX = 500;
const HEAT_MAP_BASE_HEIGHT_PX = 300;
const HEAT_MAP_BASE_CELL_PX = 20;

/**
 * Approximating the heatmap width by the number of genes.
 * This is used to make sure the table cell size stays the same regardless
 * of the number of genes selected.
 */
export function getHeatmapWidth(genes: State["selectedGenes"]): number {
  return HEAT_MAP_BASE_WIDTH_PX + HEAT_MAP_BASE_CELL_PX * genes.length;
}

/**
 * Approximating the heatmap height by the number of cells.
 */
export function getHeatmapHeight(cellTypes: CellTypeSummary[]): number {
  return HEAT_MAP_BASE_HEIGHT_PX + HEAT_MAP_BASE_CELL_PX * cellTypes.length;
}

/**
 * Value format: `${id}~${tissue}~${name}`
 * We need to encode cell type metadata here, so we can use it in onClick event
 */
export function getCellTypeNames(cellTypes: CellTypeSummary[]): string[] {
  return cellTypes.map(({ id, name }) => {
    // TODO(thuang): Add multiple tissues support.
    return `${id}~lung~${name}`;
  });
}

/**
 * Value format: `${id}~${tissue}~${name}`
 */
export function deserializeCellTypeName(cellTypeName: string): {
  id: string;
  name: string;
  tissue: string;
} {
  const [id, tissue, name] = cellTypeName.split("~");

  return {
    id,
    name,
    tissue,
  };
}

export function isCellTypeNameATissue(cellTypeName: string): boolean {
  const { name, tissue } = deserializeCellTypeName(cellTypeName);

  return name === tissue;
}

/**
 * Value format: `${name}`
 */
export function getGeneNames(genes: State["selectedGenes"]): string[] {
  return genes;
}
