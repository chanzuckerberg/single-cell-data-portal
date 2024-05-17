export interface DifferentialExpressionRow {
  name: string;
  logFoldChange: string;
  effectSize: string;
  adjustedPValue: string;
}

export interface Props {
  setIsLoading: (isLoading: boolean) => void;
}
