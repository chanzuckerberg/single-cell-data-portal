import { memo } from "react";
import ExpressedInCells from "../ExpressedInCells";
import RelativeGeneExpression from "../RelativeGeneExpression";
import { LegendWrapper } from "./style";

export default memo(function Legend(): JSX.Element {
  return (
    <LegendWrapper>
      <RelativeGeneExpression />
      <ExpressedInCells />
    </LegendWrapper>
  );
});
