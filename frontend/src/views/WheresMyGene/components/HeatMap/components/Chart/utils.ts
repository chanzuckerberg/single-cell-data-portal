import { interpolateMagma } from "d3-scale-chromatic";
import {
  CellTypeSummary,
  GeneExpressionSummary,
} from "src/views/WheresMyGene/common/types";
import { ChartFormat } from "./components/Chart/hooks/utils";
import { DefaultLabelFormatterCallbackParams, EChartsOption } from "echarts";
import { useMemo } from "react";

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
        oldRange > 0
          ? (meanExpression - scaledMeanExpressionMin) / oldRange
          : 1.0;

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
        percentage: percentage ?? tissuePercentage,
        scaledMeanExpression,
        tissuePercentage,
      };
    });
  }
}

export function useChartItemStyle(isScaled: boolean, maxExpression: number) {
  return useMemo(() => {
    return {
      color(props: DefaultLabelFormatterCallbackParams) {
        const { scaledMeanExpression, meanExpression } = props.data as {
          meanExpression: number;
          scaledMeanExpression: number;
        };

        const expressionValue = isScaled
          ? scaledMeanExpression
          : meanExpression / maxExpression;

        return interpolateMagma(1 - expressionValue);
      },
    };
  }, [isScaled, maxExpression]);
}

export function symbolSize(props: { percentage: number }) {
  const { percentage } = props;
  return convertPercentageToDiameter(percentage);
}

const Y_AXIS_BOTTOM_PADDING = "10px";

export function grid(
  defaultOption: EChartsOption["grid"]
): EChartsOption["grid"] {
  return {
    ...defaultOption,
    bottom: Y_AXIS_BOTTOM_PADDING,
  };
}

const MAX_FIRST_PART_LENGTH_PX = 16;

export function convertPercentageToDiameter(percentage: number): number {
  const maxRadius = MAX_FIRST_PART_LENGTH_PX / 2;

  const RADIUS_OFFSET = 0.2;

  const baseRadius = RADIUS_OFFSET * (MAX_FIRST_PART_LENGTH_PX - RADIUS_OFFSET);

  const radius = Math.sqrt(
    percentage * (maxRadius - RADIUS_OFFSET) ** 2 + baseRadius
  );

  return Math.round(2 * radius);
}
