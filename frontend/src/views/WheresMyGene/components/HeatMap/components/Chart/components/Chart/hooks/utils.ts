import { ViewId } from "src/views/WheresMyGene/common/types";

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
