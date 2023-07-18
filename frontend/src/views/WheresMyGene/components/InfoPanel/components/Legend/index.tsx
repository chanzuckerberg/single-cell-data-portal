import { Dispatch, memo, MouseEventHandler, SetStateAction } from "react";
import {
  FilterDimensions,
  OntologyTerm,
} from "src/common/queries/wheresMyGene";
import { CellType, ChartProps } from "src/views/WheresMyGene/common/types";
import SaveExport from "../../../GeneSearchBar/components/SaveExport";
import SaveExportV2 from "src/views/WheresMyGeneV2/components/GeneSearchBar/components/SaveExport";
import ShareButton from "../../../GeneSearchBar/components/ShareButton";
import ShareButtonV2 from "src/views/WheresMyGeneV2/components/GeneSearchBar/components/ShareButton";
import SourceDataButton from "../../../GeneSearchBar/components/SourceDataButton";
import ExpressedInCells from "../ExpressedInCells";
import RelativeGeneExpression from "../RelativeGeneExpression";
import { LegendWrapper } from "./style";
import {
  EMPTY_ARRAY,
  EMPTY_OBJECT,
  EMPTY_SET,
} from "src/common/constants/utils";

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
  selectedTissues,
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
      {selectedTissues ? (
        <>
          <SaveExport
            selectedTissues={selectedTissues}
            selectedGenes={selectedGenes}
            selectedCellTypes={selectedCellTypes}
            setDownloadStatus={setDownloadStatus}
            setEchartsRendererMode={setEchartsRendererMode}
            allChartProps={allChartProps}
            availableFilters={availableFilters}
          />
          <ShareButton />
        </>
      ) : (
        <>
          <SaveExportV2
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
          <ShareButtonV2 />
        </>
      )}
      <SourceDataButton
        handleRightSidebarButtonClick={handleRightSidebarButtonClick}
      />
      <RelativeGeneExpression isScaled={isScaled} />
      <ExpressedInCells />
    </LegendWrapper>
  );
});
