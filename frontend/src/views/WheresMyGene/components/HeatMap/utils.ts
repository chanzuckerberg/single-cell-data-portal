import { interpolateYlOrRd } from "d3-scale-chromatic";
import {
  DefaultLabelFormatterCallbackParams,
  EChartsOption,
  ScatterSeriesOption,
} from "echarts";
import { State } from "../../common/store";
import { CellTypeSummary, Tissue } from "../../common/types";
import ReplaySVG from "./icons/replay.svg";

export const MAX_FIRST_PART_LENGTH_PX = 16;
export const X_AXIS_CHART_HEIGHT = "200px";
export const Y_AXIS_CHART_WIDTH = "300px";

const REPLAY_SVG_WIDTH_PX = 12;

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
        const { name: cellTypeMetadata, data } = props as unknown as {
          name: string;
          data: ChartFormat;
        };

        if (!data) return "";

        const { geneIndex, percentage, meanExpression, scaledMeanExpression } =
          data;

        const { name } = deserializeCellTypeMetadata(
          cellTypeMetadata as CellTypeMetadata
        );

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
        axisLabel: { fontSize: 0 },
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
        axisLabel: { fontSize: 0 },
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
  tissuesWithDeletedCellTypes: string[];
}

export function createYAxisOptions({
  geneNames,
  cellTypeNames,
  commonOptions,
  cellTypeIdsToDelete,
  tissuesWithDeletedCellTypes,
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
          formatter(value) {
            const { name } = deserializeCellTypeMetadata(
              value as CellTypeMetadata
            );

            const isTissue = checkIsTissue(value as CellTypeMetadata);
            const tissueHasDeletedCellTypes =
              tissuesWithDeletedCellTypes.includes(name);

            if (isTissue) {
              return tissueHasDeletedCellTypes
                ? `{tissue|${capitalize(name)}}{reset|}`
                : `{tissue|${capitalize(name)}}`;
            }

            return cellTypeIdsToDelete.includes(value)
              ? `{selected|${name}}`
              : name;
          },
          rich: {
            reset: {
              backgroundColor: {
                image: ReplaySVG.src,
              },
              width: REPLAY_SVG_WIDTH_PX,
            },
            selected: {
              color: "red",
              fontWeight: "bold",
            },
            tissue: {
              align: "left",
              color: "black",
              fontSize: 13,
              fontWeight: "bold",
            },
          },
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
 */
export type CellTypeMetadata =
  `${CellTypeSummary["id"]}~${Tissue}~${CellTypeSummary["name"]}`;

/**
 * We need to encode cell type metadata here, so we can use it in onClick event
 */
export function getAllSerializedCellTypeMetadata(
  cellTypes: CellTypeSummary[]
): CellTypeMetadata[] {
  return cellTypes.map(({ id, name, tissue }) => {
    return `${id}~${tissue}~${name}` as CellTypeMetadata;
  });
}

export function deserializeCellTypeMetadata(
  cellTypeMetadata: CellTypeMetadata
): {
  id: string;
  name: string;
  tissue: string;
} {
  const [id, tissue, name] = cellTypeMetadata.split("~");

  return {
    id,
    name,
    tissue,
  };
}

export function checkIsTissue(cellTypeMetadata: CellTypeMetadata): boolean {
  const { name, tissue } = deserializeCellTypeMetadata(cellTypeMetadata);

  return name === tissue;
}

/**
 * Value format: `${name}`
 */
export function getGeneNames(genes: State["selectedGenes"]): string[] {
  return genes;
}

function capitalize(str: string): string {
  return str.charAt(0).toUpperCase() + str.slice(1);
}
