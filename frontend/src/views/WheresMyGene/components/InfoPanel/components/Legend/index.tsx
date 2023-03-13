import { Dispatch, memo, MouseEventHandler } from "react";
import {
  FilterDimensions,
  OntologyTerm,
} from "src/common/queries/wheresMyGene";
import { CellType } from "src/views/WheresMyGene/common/types";
import SaveImage from "../../../GeneSearchBar/components/SaveImage";
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
  setIsDownloading: (isDownloading: boolean) => void;
  setEchartsRendererMode: Dispatch<React.SetStateAction<"canvas" | "svg">>;
  allChartProps: { [tissue: string]: ChartProps };
  availableFilters: Partial<FilterDimensions>;
  availableOrganisms: OntologyTerm[];
}

export default memo(function Legend({
  isScaled,
  handleRightSidebarButtonClick,
  selectedTissues,
  selectedGenes,
  selectedCellTypes,
  setIsDownloading,
  setEchartsRendererMode,
  allChartProps,
  availableFilters,
  availableOrganisms,
}: Props): JSX.Element {
  return (
    <LegendWrapper>
      <SaveImage
        selectedTissues={selectedTissues}
        selectedGenes={selectedGenes}
        selectedCellTypes={selectedCellTypes}
        setIsDownloading={setIsDownloading}
        setEchartsRendererMode={setEchartsRendererMode}
        allChartProps={allChartProps}
        availableFilters={availableFilters}
        availableOrganisms={availableOrganisms}
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
