import React from "react";
import { useConnect } from "./connect";
import { CardWrapper, CardContent } from "./style";
import { Props } from "./types";
const InterpretationCard = ({
  isVisible,
  isQueryGroup1,
  differentialExpressionResults,
  setIsLoadingInterpret,
  setIsVisible,
}: Props) => {
  const { analysisMessage } = useConnect({
    isQueryGroup1,
    differentialExpressionResults,
    setIsLoadingInterpret,
  });

  return (
    <>
      {isVisible && (
        <CardWrapper>
          <CardContent>{analysisMessage}</CardContent>
        </CardWrapper>
      )}
      <button onClick={() => setIsVisible((prev) => !prev)}>
        {isVisible ? "Hide Analysis" : "Show Analysis"}
      </button>
    </>
  );
};

export default InterpretationCard;
