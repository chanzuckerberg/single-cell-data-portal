import { Dispatch, MouseEventHandler, SetStateAction } from "react";
import {
  FilterDimensions,
  OntologyTerm,
} from "src/common/queries/wheresMyGene";
import { CellType, ChartProps } from "src/views/WheresMyGeneV2/common/types";

export interface Props {
  isScaled: boolean;
  handleRightSidebarButtonClick: MouseEventHandler<HTMLButtonElement>;
  selectedTissues?: Array<string>;
  selectedGenes: Array<string>;
  selectedCellTypes: { [tissue: string]: CellType[] };
  setDownloadStatus: Dispatch<
    SetStateAction<{
      isLoading: boolean;
    }>
  >;
  setEchartsRendererMode: Dispatch<SetStateAction<"canvas" | "svg">>;
  allChartProps: { [tissue: string]: ChartProps };
  availableFilters: Partial<FilterDimensions>;
  tissues?: { [name: string]: OntologyTerm };
  expandedTissueIds?: string[];
  filteredCellTypes?: string[];
  maxExpression: number;
}
