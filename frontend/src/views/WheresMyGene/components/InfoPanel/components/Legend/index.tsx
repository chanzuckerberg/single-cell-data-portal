import { Dispatch, memo, MouseEventHandler, SetStateAction } from "react";
import {
  FilterDimensions,
  OntologyTerm,
} from "src/common/queries/wheresMyGene";
import { CellType, ChartProps } from "src/views/WheresMyGene/common/types";
import CitationButton from "src/views/WheresMyGeneV2/components/GeneSearchBar/components/CitationButton";
import SaveExport from "src/views/WheresMyGeneV2/components/GeneSearchBar/components/SaveExport";
import ShareButton from "src/views/WheresMyGeneV2/components/GeneSearchBar/components/ShareButton";
import ExpressedInCells from "../ExpressedInCells";
import RelativeGeneExpression from "../RelativeGeneExpression";
import { LegendWrapper } from "./style";
import {
  EMPTY_ARRAY,
  EMPTY_OBJECT,
  EMPTY_SET,
} from "src/common/constants/utils";
import SourceDataButton from "src/views/WheresMyGeneV2/components/GeneSearchBar/components/SourceDataButton";

interface Props {
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
  expandedTissues?: Set<string>;
  filteredCellTypes?: string[];
}

export default memo(function Legend({
  isScaled,
  handleRightSidebarButtonClick,
  selectedGenes,
  selectedCellTypes,
  setDownloadStatus,
  setEchartsRendererMode,
  allChartProps,
  availableFilters,
  tissues,
  expandedTissues,
  filteredCellTypes,
}: Props): JSX.Element {
  return (
    <LegendWrapper data-testid="legend-wrapper">
      <SaveExport
        selectedGenes={selectedGenes}
        selectedCellTypes={selectedCellTypes}
        setDownloadStatus={setDownloadStatus}
        setEchartsRendererMode={setEchartsRendererMode}
        allChartProps={allChartProps}
        availableFilters={availableFilters}
        tissues={tissues || EMPTY_OBJECT}
        expandedTissues={expandedTissues ?? (EMPTY_SET as Set<string>)}
        filteredCellTypes={filteredCellTypes ?? EMPTY_ARRAY}
      />
      <ShareButton />
      <CitationButton />
      <SourceDataButton
        handleRightSidebarButtonClick={handleRightSidebarButtonClick}
      />
      <RelativeGeneExpression isScaled={isScaled} />
      <ExpressedInCells />
    </LegendWrapper>
  );
});
