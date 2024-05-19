import { SetStateAction, Dispatch } from "react";
import { DifferentialExpressionResultToInterpret } from "src/common/queries/differentialExpression";

export interface Props {
  isQueryGroup1: boolean;
  differentialExpressionResults: DifferentialExpressionResultToInterpret[];
  setIsLoadingInterpret: Dispatch<SetStateAction<boolean>>;
  setIsVisible: Dispatch<SetStateAction<boolean>>;
}
