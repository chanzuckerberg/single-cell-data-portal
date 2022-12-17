import { memo, MouseEventHandler } from "react";
import { State } from "src/views/WheresMyGene/common/store";
import { CellType } from "src/views/WheresMyGene/common/types";
// import SaveImage from "../../../GeneSearchBar/components/SaveImage";
import SourceDataButton from "../../../GeneSearchBar/components/SourceDataButton";
import ExpressedInCells from "../ExpressedInCells";
import RelativeGeneExpression from "../RelativeGeneExpression";
import { LegendWrapper } from "./style";

interface Props {
  isScaled: boolean;
  handleRightSidebarButtonClick: MouseEventHandler<HTMLButtonElement>;
  selectedTissues: Array<string>;
  selectedGenes: State["selectedGenes"];
  selectedCellTypes: { [tissue: string]: CellType[] };
}

export default memo(function Legend({
  isScaled,
  handleRightSidebarButtonClick,
  // selectedTissues,
  // selectedGenes,
  // selectedCellTypes,
}: Props): JSX.Element {
  return (
    <LegendWrapper>
      {/*<SaveImage
        selectedTissues={selectedTissues}
        selectedGenes={selectedGenes}
        selectedCellTypes={selectedCellTypes}
      />*/}

      <SourceDataButton
        handleRightSidebarButtonClick={handleRightSidebarButtonClick}
      />
      <RelativeGeneExpression isScaled={isScaled} />
      <ExpressedInCells />
    </LegendWrapper>
  );
});
