import React from "react";
import { useConnect } from "./connect";
import { CardWrapper, CardContent } from "./style";
import { Props } from "./types";
const InterpretationCard = ({
  isQueryGroup1,
  differentialExpressionResults,
}: Props) => {
  const { analysisMessage, isCardVisible, toggleCardVisibility } = useConnect({
    isQueryGroup1,
    differentialExpressionResults,
  });

  return (
    <>
      {isCardVisible && (
        <CardWrapper>
          <CardContent>{analysisMessage}</CardContent>
        </CardWrapper>
      )}
      <button onClick={toggleCardVisibility}>
        {isCardVisible ? "Hide Analysis" : "Show Analysis"}
      </button>
    </>
  );
};

export default InterpretationCard;
