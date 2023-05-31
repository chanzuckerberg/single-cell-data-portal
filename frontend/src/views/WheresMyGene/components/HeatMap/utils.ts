import { interpolateMagma } from "d3-scale-chromatic";
import {
  DatasetComponentOption,
  DefaultLabelFormatterCallbackParams,
  EChartsOption,
  ScatterSeriesOption,
} from "echarts";
import { EMPTY_ARRAY } from "src/common/constants/utils";
import {
  getOntologyTermIdFromCellTypeViewId,
  getOptionIdFromCellTypeViewId,
  OptionId,
} from "src/common/queries/wheresMyGene";
import { INITIAL_STATE } from "../../common/store";
import {
  CellType,
  CellTypeSummary,
  GeneExpressionSummary,
  Tissue,
  ViewId,
} from "../../common/types";
import { TISSUE_CELL_TYPE_DIVIDER } from "./hooks/useSortedGeneNames";

export const MAX_FIRST_PART_LENGTH_PX = 16;
export const X_AXIS_HOVER_CONTAINER_HEIGHT_PX = 40;
export const X_AXIS_CHART_HEIGHT_PX = INITIAL_STATE.xAxisHeight;
export const Y_AXIS_CHART_WIDTH_PX = 300;

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

/**
 * Used to calculate text pixel widths. Should be only created once.
 */
const CTX =
  (typeof document !== "undefined" &&
    document.createElement("canvas").getContext("2d")) ||
  null;

/**
 * Formats and truncates the cell type name to a given width
 *
 * @param name The text to truncate
 * @param maxWidth The max width in pixels the string should be
 * @param font The font family and font size as a string. Ex. "bold 12px sans-serif"
 * @returns The string fixed to a certain pixel width
 */
export function formatLabel(
  name: string,
  maxWidth: number,
  font: string
): {
  text: string;
  length: number;
} {
  // failsafe, should never be falsy
  if (!CTX) return { text: name, length: 0 };

  CTX.font = font;
  const ellipsisWidth = CTX.measureText("...").width;

  const fullWidth = CTX.measureText(name).width;

  if (fullWidth <= maxWidth) {
    return {
      text: name,
      length: fullWidth,
    };
  }

  const labelHalfWidth = (maxWidth - ellipsisWidth) / 2;

  const firstHalf = getFixedWidth(name, labelHalfWidth, font);
  const secondHalf = getFixedWidth(name, labelHalfWidth, font, true);

  const formattedLabel = firstHalf + "..." + secondHalf;

  return {
    text: formattedLabel,
    length: CTX.measureText(formattedLabel).width,
  };
}

/**
 * Truncates the string to a given width
 *
 * @param text The text to truncate
 * @param maxWidth The max width in pixels the string should be
 * @param font The font family and font size as a string. Ex. "bold 12px sans-serif"
 * @param reverse Whether to truncate the end or beginning of the string
 * @returns The string fixed to a certain pixel width
 */
export function getFixedWidth(
  text: string,
  maxWidth: number,
  font: string,
  reverse = false
): string {
  // failsafe, should never be falsy
  if (!CTX) return text;

  CTX.font = font;

  if (reverse) {
    for (let i = text.length; i >= 0; i--) {
      const substring = text.substring(i - 1);
      const textWidth = CTX.measureText(substring).width;
      if (textWidth > maxWidth) {
        return text.substring(i);
      }
    }
  } else {
    for (let i = 0; i < text.length; i++) {
      const substring = text.substring(0, i + 1);
      const textWidth = CTX.measureText(substring).width;
      if (textWidth > maxWidth) {
        return text.substring(0, i);
      }
    }
  }

  return text;
}

interface CreateChartOptionsProps {
  cellTypeMetadata: string[];
  chartData: ChartFormat[];
  geneNames: string[];
  isScaled: boolean;
  heatmapWidth: number;
  heatmapHeight: number;
}

export function createChartOptions({
  cellTypeMetadata,
  chartData,
  geneNames,
  isScaled,
  heatmapWidth,
  heatmapHeight,
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
      height: `${heatmapHeight}px`,
      left: "0px",
      top: "0px",
      // (atarashansky): this is the key change to align x and y axis
      // labels to fixed spacings
      width: `${heatmapWidth}px`,
    },
    series: [
      {
        ...COMMON_SERIES,
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
export interface ChartFormat {
  cellTypeIndex: number;
  geneIndex: number;
  percentage: number;
  meanExpression: number;
  scaledMeanExpression: number;
  tissuePercentage: number;
  expressedCellCount: number;
  id: `${ViewId}-${string}`;
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
        (cellTypeSummary) => cellTypeSummary.viewId === dataPoint.viewId
      );

      const id = `${dataPoint.viewId}-${geneName}` as ChartFormat["id"];

      return {
        cellTypeIndex,
        expressedCellCount,
        geneIndex,
        id,
        meanExpression,
        percentage,
        scaledMeanExpression,
        tissuePercentage,
      };
    });
  }
}

export const HEAT_MAP_BASE_HEIGHT_PX = 300;
export const HEAT_MAP_BASE_CELL_PX = 20;
export const HEAT_MAP_BASE_CELL_WIDTH_PX = 20;

/**
 * Approximating the heatmap width by the number of genes.
 * This is used to make sure the table cell size stays the same regardless
 * of the number of genes selected.
 */
export function getHeatmapWidth(
  genes: (GeneExpressionSummary | undefined)[] | string[] = EMPTY_ARRAY
): number {
  return HEAT_MAP_BASE_CELL_WIDTH_PX * genes.length;
}

/**
 * Approximating the heatmap height by the number of cells.
 */
export function getHeatmapHeight(cellTypes: CellType[] = EMPTY_ARRAY): number {
  return HEAT_MAP_BASE_CELL_PX * cellTypes.length;
}

/**
 * Value format: `${viewId}~${tissue}~${name}~${order}`
 */
export type CellTypeMetadata =
  `${CellTypeSummary["viewId"]}~${Tissue}~${CellTypeSummary["name"]}~${number}~${number}`;

/**
 * We need to encode cell type metadata here, so we can use it in onClick event
 */
export function getAllSerializedCellTypeMetadata(
  cellTypes: CellType[],
  tissue: Tissue
): CellTypeMetadata[] {
  return cellTypes.map(({ viewId, name, order, total_count }) => {
    return `${viewId}~${tissue}~${name}~${order}~${total_count}` as CellTypeMetadata;
  });
}

export function deserializeCellTypeMetadata(
  cellTypeMetadata: CellTypeMetadata
): {
  viewId: CellType["viewId"];
  name: CellType["name"];
  tissue: Tissue;
  order: number;
  total_count: number;
  id: string;
  optionId: OptionId;
} {
  const [rawViewId, tissue, name, order, total_count] = cellTypeMetadata.split(
    TISSUE_CELL_TYPE_DIVIDER
  );

  // (thuang): Typescript limitation that split() returns string[]
  const viewId = rawViewId as CellType["viewId"];

  const id = getOntologyTermIdFromCellTypeViewId(viewId);
  const optionId = getOptionIdFromCellTypeViewId(viewId);

  return {
    id,
    name,
    optionId,
    order: Number(order),
    tissue,
    total_count: Number(total_count),
    viewId,
  };
}

/**
 * Value format: `${name}`
 */
export function getGeneNames(
  genes: ({ name: string } | undefined)[]
): string[] {
  return genes.map((gene) => gene?.name || "");
}
