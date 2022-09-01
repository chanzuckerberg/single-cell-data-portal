import { memo } from "react";
import ExpressedInCells from "../ExpressedInCells";
import RelativeGeneExpression from "../RelativeGeneExpression";
import { LegendWrapper } from "./style";

interface Props {
  handleIsScaledChange: () => void;
  isScaled: boolean;
  showScaled?: boolean;
}

export default memo(function InfoPanel({
  handleIsScaledChange,
  isScaled,
  showScaled,
}: Props): JSX.Element {
  return (
    <LegendWrapper>
      <RelativeGeneExpression
        isScaled={isScaled}
        handleIsScaledChange={handleIsScaledChange}
        showScaled={showScaled}
      />
      <ExpressedInCells />
    </LegendWrapper>
  );
});
