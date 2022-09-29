import { memo, MouseEventHandler } from "react";
import SaveImage from "../../../GeneSearchBar/components/SaveImage";
import SourceDataButton from "../../../GeneSearchBar/components/SourceDataButton";
import ExpressedInCells from "../ExpressedInCells";
import RelativeGeneExpression from "../RelativeGeneExpression";
import { LegendWrapper } from "./style";

interface Props {
  isScaled: boolean;
  handleRightSidebarButtonClick: MouseEventHandler<HTMLButtonElement>;
  selectedTissues: Array<string>;
  selectedGenes: Array<string>;
}

export default memo(function Legend({
  isScaled,
  handleRightSidebarButtonClick,
  selectedTissues,
  selectedGenes,
}: Props): JSX.Element {
  return (
    <LegendWrapper>
      <SaveImage
        selectedTissues={selectedTissues}
        selectedGenes={selectedGenes}
      />
      <SourceDataButton
        handleRightSidebarButtonClick={handleRightSidebarButtonClick}
      />
      <RelativeGeneExpression isScaled={isScaled} />
      <ExpressedInCells />
    </LegendWrapper>
  );
});
