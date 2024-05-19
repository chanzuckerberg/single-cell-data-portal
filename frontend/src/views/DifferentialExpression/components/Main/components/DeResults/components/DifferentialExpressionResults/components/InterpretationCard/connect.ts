import { useState, useEffect } from "react";
import { useInterpretDeResults } from "src/common/queries/differentialExpression";
import { Props } from "./types";

export const useConnect = ({
  isQueryGroup1,
  differentialExpressionResults,
  setIsLoadingInterpret,
}: Props) => {
  const [isCardVisible, setIsCardVisible] = useState(false);
  const {
    data: { message: analysisMessage },
    isLoading,
  } = useInterpretDeResults({
    isQueryGroup1,
    differentialExpressionResults,
  });

  useEffect(() => {
    setIsLoadingInterpret(isLoading);
  }, [isLoading, setIsLoadingInterpret]);

  const toggleCardVisibility = () => {
    setIsCardVisible((prev) => !prev);
  };

  return {
    analysisMessage,
    isCardVisible,
    toggleCardVisibility,
  };
};
