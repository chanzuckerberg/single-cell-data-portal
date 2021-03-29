export interface RowInfo {
  row: string[];
  rowIndex: number;
  errorMessages: Set<string>;
  lastGeneSetName: string;
  result: Result;
  headerRow: string[] | null;
}

export type Result = Map<string, GeneSet>;

export const COMMENT_SYMBOL = "#";

export interface GeneSet {
  gene_set_name: string;
  gene_set_description: string;
  genes: Map<string, Gene>;
}

export interface Gene {
  gene_symbol: string;
  gene_description?: string;
  additional_params?: {
    [key: string]: string;
  };
}

export interface OutputGeneSet extends Omit<GeneSet, "genes"> {
  genes: Gene[];
}

export type Output = OutputGeneSet[];
