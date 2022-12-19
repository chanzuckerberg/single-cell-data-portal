import { memo, MouseEventHandler } from "react";
import { State } from "src/views/WheresMyGene/common/store";
import { CellType } from "src/views/WheresMyGene/common/types";
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
}: 
Props): JSX.Element {
  return (
    <LegendWrapper>
      <SourceDataButton
        handleRightSidebarButtonClick={handleRightSidebarButtonClick}
      />
      <RelativeGeneExpression isScaled={isScaled} />
      <ExpressedInCells />
    </LegendWrapper>
  );
});
