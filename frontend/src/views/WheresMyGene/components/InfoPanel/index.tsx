import { memo } from "react";
import Methodology from "./components/Methodology";
import SourceData from "./components/SourceData";

interface Props {
  isScaled: boolean;
}

export default memo(function InfoPanel(): JSX.Element {
  return (
    <>
      <SourceData />
      <Methodology />
    </>
  );
});
