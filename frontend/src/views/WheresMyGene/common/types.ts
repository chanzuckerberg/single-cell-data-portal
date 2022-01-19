import { CSSProperties } from "react";

/** This object holds a cell type and its genes. */
export interface CellTypeAndGenes {
  id: string;
  name: string;
  expressions?: {
    [geneName: string]: GeneExpression;
  };
}

export interface CellType {
  id: string;
  name: string;
}

/**
 * This object holds a gene and all the cell types that has this
 * genetic expression.
 */
export interface RawGeneExpression {
  gene_name: string;
  cell_types: GeneExpression[];
}

/**
 * This object describes the cell type id and the cell's genetic
 * expression metadata
 */
export interface GeneExpression {
  /** cellTypeId */
  id: string;
  /** percentage of the cells that have this genetic expression */
  pc: number;
  /** mean expression of the cells that have this genetic expression */
  me: number;
}

export interface Gene {
  id: string;
  name: string;
  style?: CSSProperties;
  cell_types?: GeneExpression[];
}
