import { interpolateMagma } from "d3-scale-chromatic";
import {
  DatasetComponentOption,
  DefaultLabelFormatterCallbackParams,
  EChartsOption,
  ScatterSeriesOption,
} from "echarts";
import { EMPTY_ARRAY } from "src/common/constants/utils";
import { LIGHT_GRAY } from "src/components/common/theme";
import { State } from "../../common/store";
import {
  CellType,
  CellTypeSummary,
  GeneExpressionSummary,
  Tissue,
} from "../../common/types";

export const MAX_FIRST_PART_LENGTH_PX = 16;
export const X_AXIS_CHART_HEIGHT_PX = 80;
export const Y_AXIS_CHART_WIDTH_PX = 300;

const MAX_DEPTH = 2;

const Y_AXIS_BOTTOM_PADDING = "10px";

const COMMON_OPTIONS = {
  animation: false,
  hoverLayerThreshold: 10,
  progressive: 1e6,
};

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

const MAX_MEAN_EXPRESSION_VALUE = 6;

interface CreateChartOptionsProps {
  cellTypeMetadata: string[];
  chartData: ChartFormat[];
  geneNames: string[];
  isScaled: boolean;
}

export function createChartOptions({
  cellTypeMetadata,
  chartData,
  geneNames,
  isScaled,
}: CreateChartOptionsProps): EChartsOption {
  return {
    ...COMMON_OPTIONS,
    axisPointer: {
      label: { show: false },
      link: [
        {
          xAxisIndex: "all",
        },
      ],
      show: true,
      triggerOn: "mousemove",
    },
    dataset: {
      source: chartData as DatasetComponentOption["source"],
    },
    grid: {
      bottom: Y_AXIS_BOTTOM_PADDING,
    },
    series: [
      {
        ...COMMON_SERIES,
        symbolOffset: [0,-60],
        itemStyle: {
          color(props: DefaultLabelFormatterCallbackParams) {
            const { scaledMeanExpression, meanExpression } = props.data as {
              meanExpression: number;
              scaledMeanExpression: number;
            };

            const expressionValue = isScaled
              ? scaledMeanExpression
              : meanExpression / MAX_MEAN_EXPRESSION_VALUE;

            return interpolateMagma(1 - expressionValue);
          },
        },
        symbolSize: function (props: { percentage: number }) {
          const { percentage } = props;

          return convertPercentageToDiameter(percentage);
        },
      },
    ],
    xAxis: [
      {
        axisLabel: { fontSize: 0 },
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
        axisLabel: { fontSize: 0 },
        axisLine: {
          show: false,
        },
        axisTick: {
          show: false,
        },
        boundaryGap: true,
        data: cellTypeMetadata,
        splitLine: {
          show: false,
        },
      },
    ],
  };
}

export function convertPercentageToDiameter(percentage: number): number {
  const maxRadius = MAX_FIRST_PART_LENGTH_PX / 2;

  const RADIUS_OFFSET = 0.2;

  const baseRadius = RADIUS_OFFSET * (MAX_FIRST_PART_LENGTH_PX - RADIUS_OFFSET);

  const radius = Math.sqrt(
    percentage * (maxRadius - RADIUS_OFFSET) ** 2 + baseRadius
  );

  return Math.round(2 * radius);
}

const SELECTED_STYLE = {
  backgroundColor: LIGHT_GRAY.D,
  fontWeight: "bold" as never,
  padding: 4,
};

interface CreateXAxisOptionsProps {
  geneNames: string[];
  genesToDelete: string[];
  noSelect?: boolean;
}

export function createXAxisOptions({
  geneNames,
  genesToDelete,
  noSelect
}: CreateXAxisOptionsProps): EChartsOption {
  return {
    ...COMMON_OPTIONS,
    grid: {
      bottom: "0",
      top: "240px",
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
            return genesToDelete.includes(value) && !noSelect
              ? `{selected|${value}}`
              : value;
          },
          rich: {
            selected: SELECTED_STYLE,
          },
          rotate: 270,
          verticalAlign: "middle",
          width: 200,
        },
        axisTick: {
          show: false,
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
        //name: "Cell Count",
        //nameLocation: "end",
        //nameRotate: 270,
        //offset: 10,
        //nameGap: 8,
        axisLabel: { fontSize: 0 },
        axisLine: {
          show: false,
        },
        axisTick: {
          show: false,
        },
        splitLine: {
          show: false,
        },
      },
    ],
  };
}

interface CreateYAxisOptionsProps {
  cellTypeMetadata: CellTypeMetadata[];
  cellTypeIdsToDelete: string[];
}

export function createYAxisOptions({
  cellTypeMetadata,
  cellTypeIdsToDelete,
}: CreateYAxisOptionsProps): EChartsOption {
  return {
    ...COMMON_OPTIONS,
    grid: {
      bottom: Y_AXIS_BOTTOM_PADDING,
      left: "20px",
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
        splitLine: {
          show: false,
        },
        type: "category",
      },
    ],
    yAxis: [
      {
        axisLabel: {
          align: "left",
          formatter(value: number | string) {
            const { name, depth = 0 } = deserializeCellTypeMetadata(
              value as CellTypeMetadata
            );

            const displayDepth = Math.min(depth, MAX_DEPTH);

            const paddedName = " ".repeat(displayDepth * 8) + name;

            return cellTypeIdsToDelete.includes(value as string)
              ? `{selected|${paddedName}}`
              : paddedName;
          },
          // Turn off type checking here, because ecahrts' type is wrong
          ["overflow" as string]: "truncate",
          rich: {
            selected: SELECTED_STYLE,
          },
          width: 230,
        },
        axisLine: {
          show: false,
        },
        axisTick: {
          show: false,
        },
        boundaryGap: true,
        data: cellTypeMetadata,
        splitLine: {
          show: false,
        },
        triggerEvent: true,
      },
      {
        axisLabel: {
          formatter(value: number | string) {
            const { total_count } = deserializeCellTypeMetadata(
              value as CellTypeMetadata
            );
            return total_count > 10000 ? ">10k" : `${total_count}`;
          },

          rich: {
            selected: SELECTED_STYLE,
          },
          align: "right",
        },
        position: "right",
        offset: -10,
        axisLine: {
          show: false,
        },
        axisTick: {
          show: false,
        },
        boundaryGap: true,
        data: cellTypeMetadata,
        splitLine: {
          show: false,
        },
      },
    ],
  };
}

export interface ChartFormat {
  cellTypeIndex: number;
  geneIndex: number;
  percentage: number;
  id: string;
  meanExpression: number;
  scaledMeanExpression: number;
  tissuePercentage: number;
  expressedCellCount: number;
}

export function dataToChartFormat({
  cellTypeSummaries,
  genes,
  scaledMeanExpressionMax,
  scaledMeanExpressionMin,
}: {
  cellTypeSummaries: CellTypeSummary[];
  genes: (GeneExpressionSummary | undefined)[];
  scaledMeanExpressionMax: number;
  scaledMeanExpressionMin: number;
}): ChartFormat[] {
  const oldRange = scaledMeanExpressionMax - scaledMeanExpressionMin;

  const result = cellTypeSummaries.flatMap((dataPoint) => {
    return toChartFormat(dataPoint);
  });

  return result;

  function toChartFormat(dataPoint: CellTypeSummary): ChartFormat[] {
    const { geneExpressions } = dataPoint;

    if (!geneExpressions) return [];

    return Object.entries(geneExpressions).map(([geneName, geneExpression]) => {
      const {
        percentage,
        meanExpression,
        tissuePercentage,
        expressedCellCount,
      } = geneExpression;

      const scaledMeanExpression =
        (meanExpression - scaledMeanExpressionMin) / oldRange;

      const geneIndex = genes.findIndex((gene) => gene?.name === geneName);

      const cellTypeIndex = cellTypeSummaries.findIndex(
        (cellTypeSummary) => cellTypeSummary.id === dataPoint.id
      );

      const id = `${dataPoint.id}-${geneName}`;

      return {
        cellTypeIndex,
        expressedCellCount,
        geneIndex,
        id,
        meanExpression,
        percentage,
        scaledMeanExpression,
        tissuePercentage,
      } as ChartFormat;
    });
  }
}

export const HEAT_MAP_BASE_HEIGHT_PX = 300;
export const HEAT_MAP_BASE_CELL_PX = 20;
export const HEAT_MAP_BASE_CELL_WIDTH_PX = 50;

/**
 * Approximating the heatmap width by the number of genes.
 * This is used to make sure the table cell size stays the same regardless
 * of the number of genes selected.
 */
export function getHeatmapWidth(
  genes:
    | (GeneExpressionSummary | undefined)[]
    | State["selectedGenes"] = EMPTY_ARRAY
): number {
  return HEAT_MAP_BASE_CELL_WIDTH_PX*genes.length;
}

/**
 * Approximating the heatmap height by the number of cells.
 */
export function getHeatmapHeight(cellTypes: CellType[] = EMPTY_ARRAY): number {
  return HEAT_MAP_BASE_HEIGHT_PX + HEAT_MAP_BASE_CELL_PX * cellTypes.length;
}

/**
 * Value format: `${id}~${tissue}~${name}~${depth}`
 */
export type CellTypeMetadata =
  `${CellTypeSummary["id"]}~${Tissue}~${CellTypeSummary["name"]}~${number}~${number}`;

/**
 * We need to encode cell type metadata here, so we can use it in onClick event
 */
export function getAllSerializedCellTypeMetadata(
  cellTypes: CellType[],
  tissue: Tissue
): CellTypeMetadata[] {
  return cellTypes.map(({ id, name, depth, total_count }) => {
    return `${id}~${tissue}~${name}~${depth}~${total_count}` as CellTypeMetadata;
  });
}

export function deserializeCellTypeMetadata(
  cellTypeMetadata: CellTypeMetadata
): {
  id: string;
  name: string;
  tissue: Tissue;
  depth: number;
  total_count: number;
} {
  const [id, tissue, name, depth, total_count] = cellTypeMetadata.split("~");

  return {
    depth: Number(depth),
    id,
    name,
    tissue,
    total_count: Number(total_count),
  };
}

export function checkIsTissue(cellTypeMetadata: CellTypeMetadata): boolean {
  const { name, tissue } = deserializeCellTypeMetadata(cellTypeMetadata);

  return name === tissue;
}

/**
 * Value format: `${name}`
 */
export function getGeneNames(
  genes: ({ name: string } | undefined)[]
): string[] {
  return genes.map((gene) => gene?.name || "");
}
