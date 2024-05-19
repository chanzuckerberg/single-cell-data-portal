import { useEffect } from "react";
import { useInterpretDeResults } from "src/common/queries/differentialExpression";
import { Props } from "./types";

export const useConnect = ({
  isQueryGroup1,
  differentialExpressionResults,
  setIsLoadingInterpret,
}: {
  isQueryGroup1: Props["isQueryGroup1"];
  differentialExpressionResults: Props["differentialExpressionResults"];
  setIsLoadingInterpret: Props["setIsLoadingInterpret"];
}) => {
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

  return {
    analysisMessage,
    isLoading,
  };
};
