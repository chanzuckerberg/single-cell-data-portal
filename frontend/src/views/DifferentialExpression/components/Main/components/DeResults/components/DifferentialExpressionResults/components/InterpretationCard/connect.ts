import { useEffect, useMemo } from "react";
import { useInterpretDeResults } from "src/common/queries/differentialExpression";
import { Props } from "./types";
import {
  MAX_NUM_TOP_GENES_TO_INTERPRET,
  MAX_P_VALUE_TO_INTERPRET,
  MIN_EFEFCT_SIZE_TO_INTERPRET,
} from "./constants";

export const useConnect = ({
  differentialExpressionResults,
  setIsLoadingInterpret,
}: {
  differentialExpressionResults: Props["differentialExpressionResults"];
  setIsLoadingInterpret: Props["setIsLoadingInterpret"];
}) => {
  const deGenes1 = useMemo(
    () =>
      differentialExpressionResults
        .filter(
          (row) =>
            parseFloat(row.effectSize) > MIN_EFEFCT_SIZE_TO_INTERPRET &&
            parseFloat(row.adjustedPValue) < MAX_P_VALUE_TO_INTERPRET
        )
        .slice(0, MAX_NUM_TOP_GENES_TO_INTERPRET)
        .map((row) => row.name),
    [differentialExpressionResults]
  );
  const deGenes2 = useMemo(
    () =>
      differentialExpressionResults

        .filter(
          (row) =>
            parseFloat(row.effectSize) < -MIN_EFEFCT_SIZE_TO_INTERPRET &&
            parseFloat(row.adjustedPValue) < MAX_P_VALUE_TO_INTERPRET
        )
        .slice(-MAX_NUM_TOP_GENES_TO_INTERPRET)
        .reverse()
        .map((row) => row.name),
    [differentialExpressionResults]
  );
  const {
    data: { message: analysisMessage, prompt: analysisPrompt },
    isLoading,
  } = useInterpretDeResults({
    deGenes1,
    deGenes2,
  });

  useEffect(() => {
    setIsLoadingInterpret(isLoading);
  }, [isLoading, setIsLoadingInterpret]);

  return {
    analysisMessage,
    analysisPrompt,
    isLoading,
  };
};
