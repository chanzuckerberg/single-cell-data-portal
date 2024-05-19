import { useState } from "react";
import { useInterpretDeResults } from "src/common/queries/differentialExpression";
import { Props } from "./types";

export const useConnect = ({
  isQueryGroup1,
  differentialExpressionResults,
}: Props) => {
  const [isCardVisible, setIsCardVisible] = useState(false);
  const {
    data: { message: analysisMessage },
  } = useInterpretDeResults({
    isQueryGroup1,
    differentialExpressionResults,
  });

  const toggleCardVisibility = () => {
    setIsCardVisible((prev) => !prev);
  };

  return {
    analysisMessage,
    isCardVisible,
    toggleCardVisibility,
  };
};
