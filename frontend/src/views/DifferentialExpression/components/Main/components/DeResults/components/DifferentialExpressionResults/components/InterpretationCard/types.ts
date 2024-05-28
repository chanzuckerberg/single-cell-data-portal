import { SetStateAction, Dispatch } from "react";
import { DifferentialExpressionRow } from "../../../../types";

export interface Props {
  differentialExpressionResults: DifferentialExpressionRow[];
  setIsLoadingInterpret: Dispatch<SetStateAction<boolean>>;
  setIsVisible: Dispatch<SetStateAction<boolean>>;
}
