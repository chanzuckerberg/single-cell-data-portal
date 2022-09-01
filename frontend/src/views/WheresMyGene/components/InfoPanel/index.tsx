import { memo } from "react";
import ExpressedInCells from "./components/ExpressedInCells";
import Methodology from "./components/Methodology";
import RelativeGeneExpression from "./components/RelativeGeneExpression";
import SourceData from "./components/SourceData";

interface Props {
  handleIsScaledChange: () => void;
  isScaled: boolean;
}

export default memo(function InfoPanel({
  handleIsScaledChange,
  isScaled,
}: Props): JSX.Element {
  return (
    <>
      <RelativeGeneExpression
        isScaled={isScaled}
        handleIsScaledChange={handleIsScaledChange}
      />
      <ExpressedInCells />
      <Methodology />
      <SourceData />
    </>
  );
});
