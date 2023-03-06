import { Dispatch, memo, MouseEventHandler } from "react";
import { CellType } from "src/views/WheresMyGene/common/types";
import SaveImage from "../../../GeneSearchBar/components/SaveImage";
import ShareButton from "../../../GeneSearchBar/components/ShareButton";
import SourceDataButton from "../../../GeneSearchBar/components/SourceDataButton";
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
}

export default memo(function Legend({
  isScaled,
  handleRightSidebarButtonClick,
  selectedTissues,
  selectedGenes,
  selectedCellTypes,
  setIsDownloading,
  setEchartsRendererMode,
}: Props): JSX.Element {
  return (
    <LegendWrapper>
      <SaveImage
        selectedTissues={selectedTissues}
        selectedGenes={selectedGenes}
        selectedCellTypes={selectedCellTypes}
        setIsDownloading={setIsDownloading}
        setEchartsRendererMode={setEchartsRendererMode}
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
