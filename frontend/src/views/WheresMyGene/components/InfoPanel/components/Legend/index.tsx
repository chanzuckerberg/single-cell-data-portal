import { Dispatch, memo, MouseEventHandler, SetStateAction } from "react";
import { FilterDimensions } from "src/common/queries/wheresMyGene";
import { CellType } from "src/views/WheresMyGene/common/types";
import SaveExport from "../../../GeneSearchBar/components/SaveExport";
import ShareButton from "../../../GeneSearchBar/components/ShareButton";
import SourceDataButton from "../../../GeneSearchBar/components/SourceDataButton";
import { ChartProps } from "../../../HeatMap/hooks/common/types";
import ExpressedInCells from "../ExpressedInCells";
import RelativeGeneExpression from "../RelativeGeneExpression";
import { LegendWrapper } from "./style";

interface Props {
  isScaled: boolean;
  handleRightSidebarButtonClick: MouseEventHandler<HTMLButtonElement>;
  selectedTissues: Array<string>;
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
}: Props): JSX.Element {
  return (
    <LegendWrapper>
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
      <SourceDataButton
        handleRightSidebarButtonClick={handleRightSidebarButtonClick}
      />
      <RelativeGeneExpression isScaled={isScaled} />
      <ExpressedInCells />
    </LegendWrapper>
  );
});
