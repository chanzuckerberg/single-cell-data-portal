import { CSSProperties } from "react";

export interface CellTypeAndGenes {
  id: string;
  name: string;
  expressions: {
    [geneName: string]: GeneExpression;
  };
}

export interface CellType {
  id: string;
  name: string;
}

export interface RawGeneExpression {
  gene_name: string;
  cell_types: GeneExpression[];
}

export interface GeneExpression {
  id: string;
  pc: number;
  me: number;
}

export interface Gene {
  id: string;
  name: string;
  style?: CSSProperties;
  cell_types?: GeneExpression[];
}
