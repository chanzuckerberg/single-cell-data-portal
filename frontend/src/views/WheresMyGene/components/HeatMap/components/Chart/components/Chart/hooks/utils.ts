import {
  EChartsOption,
  DatasetComponentOption,
  ScatterSeriesOption,
} from "echarts";
import { ViewId } from "src/views/WheresMyGene/common/types";

export interface CreateChartOptionsProps {
  /**
   * The data array to be visualized
   * The data point object shape can be whatever you like, but it must be consistent with the `encode` option
   * For example, if the data point shape is:
   * {
   *   geneIndex: 0,
   *   cellTypeIndex: 0,
   *   percentage: 0.5
   * }
   * and you want geneIndex to be encoded to x axis and cellTypeIndex to be encoded to y axis, then make sure your encode option is:
   * encode: {
   *   x: 'geneIndex',
   *   y: 'cellTypeIndex'
   * }
   */
  data: DatasetComponentOption["source"];
  /**
   * The data for the x axis
   * For example:
   * [{ value: "gene1", textStyle: { color: "red" } }, "gene2", "gene3"]
   */
  xAxisData: CategoryAxisData;
  /**
   * The data for the y axis
   * For example:
   * [{ value: "cellType1", textStyle: { color: "red" } }, "cellType2", "cellType3"]
   */
  yAxisData: CategoryAxisData;
  width: number;
  height: number;
  /**
   * Provide a mapping of data key to x/y axis encoding
   * For example, if the data is:
   * {
   *   geneIndex: 0,
   *   cellTypeIndex: 0,
   *   percentage: 0.5
   * }
   * and we want to encode `geneIndex` to x axis and `cellTypeIndex` to y axis, then
   * encode: {
   *  x: 'geneIndex',
   *  y: 'cellTypeIndex'
   * }
   * https://echarts.apache.org/en/option.html#series-scatter.encode
   */
  encode?: {
    x: string;
    y: string;
  };
  /**
   * Customize the style of each cell item, such as color, border, opacity, etc.
   * https://echarts.apache.org/en/option.html#series-scatter.itemStyle
   */
  itemStyle?: ScatterSeriesOption["itemStyle"];
  /**
   * `symbolSize` can be set to single numbers like 10, or use an array to represent width and height. For example, [20, 10] means symbol width is 20, and height is 10.
   *
   * If size of symbols needs to be different, you can set with callback function in the following format:
   *
   * (value: Array|number, params: Object) => number|Array
   *
   * The first parameter value is the value in data, and the second parameter params is the rest parameters of data item.
   * https://echarts.apache.org/en/option.html#series-scatter.symbolSize
   */
  symbolSize?: ScatterSeriesOption["symbolSize"];
  /**
   * https://echarts.apache.org/en/option.html#grid
   */
  grid?:
    | EChartsOption["grid"]
    | ((defaultOption: EChartsOption["grid"]) => EChartsOption["grid"]);
}

export function createChartOptions({
  data,
  xAxisData,
  yAxisData,
  width,
  height,
  encode,
  itemStyle,
  symbolSize,
  grid: gridProp,
}: CreateChartOptionsProps): EChartsOption {
  const defaultGrid = {
    height: `${height}px`,
    left: 0,
    top: 0,
    // (atarashansky): this is the key change to align x and y axis
    // labels to fixed spacings
    width: `${width}px`,
  };

  const customGrid =
    typeof gridProp === "function" ? gridProp(defaultGrid) : gridProp;

  return {
    animation: false,
    /**
     * Display reference line and axis value under mouse pointer
     * https://echarts.apache.org/en/option.html#axisPointer
     */
    axisPointer: {
      label: { show: false },
      show: true,
      triggerOn: "mousemove",
    },
    dataset: {
      source: data as DatasetComponentOption["source"],
    },
    grid: customGrid || defaultGrid,
    series: [
      {
        emphasis: {
          itemStyle: {
            color: "inherit",
          },
          scale: false,
        },
        encode,
        itemStyle,
        legendHoverLink: false,
        symbolSize,
        type: "scatter",
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
        data: xAxisData,
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
        data: yAxisData,
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
  meanExpression: number;
  scaledMeanExpression: number;
  tissuePercentage: number;
  expressedCellCount: number;
  id: `${ViewId}-${string}`;
}

type OrdinalRawValue = string | number;

/**
 * (thuang): This copies echarts' CategoryAxisBaseOption["data"] type, since it's not exported
 */
type CategoryAxisData =
  | (
      | OrdinalRawValue
      | {
          value: OrdinalRawValue;
          /**
           * (thuang): This should be echarts `TextCommonOption` type, but it's not exported
           */
          textStyle?: never;
        }
    )[];
