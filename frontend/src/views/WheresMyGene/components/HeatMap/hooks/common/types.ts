import { ChartFormat } from "../../utils";

export interface ChartProps {
  chartData: ChartFormat[];
  geneNames: string[];
  cellTypeNames: string[];
}
