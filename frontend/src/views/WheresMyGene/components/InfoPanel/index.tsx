import { memo } from "react";
import ExpressedInCells from "./components/ExpressedInCells";
import Methodology from "./components/Methodology";
import RelativeGeneExpression from "./components/RelativeGeneExpression";
import SourceData from "./components/SourceData";

export default memo(function InfoPanel(): JSX.Element {
  return (
    <>
      <RelativeGeneExpression />
      <ExpressedInCells />
      <Methodology />
      <SourceData />
    </>
  );
});
