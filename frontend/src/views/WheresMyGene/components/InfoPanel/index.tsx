import { memo } from "react";
import ExpressedInCells from "./components/ExpressedInCells";
import Methodology from "./components/Methodology";
import RelativeGeneExpression from "./components/RelativeGeneExpression";
import SourceData from "./components/SourceData";

interface Props {
  isScaled: boolean;
}

export default memo(function InfoPanel({ isScaled }: Props): JSX.Element {
  return (
    <>
      <RelativeGeneExpression isScaled={isScaled} />
      <ExpressedInCells />
      <Methodology />
      <SourceData />
    </>
  );
});
