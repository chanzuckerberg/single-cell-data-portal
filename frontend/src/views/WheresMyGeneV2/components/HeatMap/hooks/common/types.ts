import { CellTypeMetadata, ChartFormat } from "../../utils";

export interface ChartProps {
  chartData: ChartFormat[];
  geneNames: string[];
  cellTypeMetadata: CellTypeMetadata[];
}
