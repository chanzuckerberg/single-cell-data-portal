import { memo } from "react";

import CitationButton from "src/views/WheresMyGeneV2/components/GeneSearchBar/components/CitationButton";
import SaveExport from "src/views/WheresMyGeneV2/components/GeneSearchBar/components/SaveExport";
import ShareButton from "src/views/WheresMyGeneV2/components/GeneSearchBar/components/ShareButton";
import ExpressedInCells from "../ExpressedInCells";
import RelativeGeneExpression from "../RelativeGeneExpression";
import { ActionsWrapper, ColorLegendWrapper, LegendWrapper } from "./style";
import { EMPTY_ARRAY, EMPTY_OBJECT } from "src/common/constants/utils";
import SourceDataButton from "src/views/WheresMyGeneV2/components/GeneSearchBar/components/SourceDataButton";
import { useConnect } from "./connect";
import { Props } from "./types";

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
  expandedTissueIds,
  filteredCellTypes,
  maxExpression,
}: Props): JSX.Element {
  const { referenceCount } = useConnect();
  return (
    <LegendWrapper data-testid="legend-wrapper">
      <ActionsWrapper>
        <SaveExport
          selectedGenes={selectedGenes}
          selectedCellTypes={selectedCellTypes}
          setDownloadStatus={setDownloadStatus}
          setEchartsRendererMode={setEchartsRendererMode}
          allChartProps={allChartProps}
          availableFilters={availableFilters}
          tissues={tissues || EMPTY_OBJECT}
          expandedTissueIds={expandedTissueIds ?? EMPTY_ARRAY}
          filteredCellTypes={filteredCellTypes ?? EMPTY_ARRAY}
        />
        <ShareButton />
        <CitationButton />
        <SourceDataButton
          handleRightSidebarButtonClick={handleRightSidebarButtonClick}
          referenceCount={referenceCount}
        />
      </ActionsWrapper>
      <ColorLegendWrapper>
        <RelativeGeneExpression
          isScaled={isScaled}
          maxExpression={maxExpression}
        />
        <ExpressedInCells />
      </ColorLegendWrapper>
    </LegendWrapper>
  );
});
