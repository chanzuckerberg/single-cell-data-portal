import { DefaultMenuSelectOption } from "czifui";
import { CSSProperties } from "react";

export interface Organism {
  id: string;
  name: string;
}

/** tissue name */
export type Tissue = string;

/** This object holds a cell type and its gene expressions. */
export interface CellTypeSummary {
  /** cellType id */
  id: string;
  /** cellType name */
  name: string;
  geneExpressions?: {
    [geneName: string]: CellTypeGeneExpressionSummaryData;
  };
}

export interface CellType {
  id: string;
  name: string;
}

/**
 * This object holds a gene and all the cell types that express this gene.
 */
export interface GeneExpressionSummary {
  /** gene name */
  name: string;
  cellTypeGeneExpressionSummaries: CellTypeGeneExpressionSummaryData[];
}

/** This is the original data shape in API response for `CellTypeGeneExpressionSummaryData` */
export interface RawCellTypeGeneExpressionSummaryData {
  /** cellTypeId */
  id: string;
  /** percentage of the current subset of cells that express this gene */
  pc: number;
  /** mean expression of the current subset of cells that express this gene */
  me: number;
}

/**
 * This object describes the cell type id and its gene expression metadata given
 * the current subset of cells.
 */
export interface CellTypeGeneExpressionSummaryData {
  /** cellTypeId */
  id: string;
  /** percentage of the current subset of cells that express this gene */
  percentage: number;
  /** mean expression of the current subset of cells that express this gene */
  meanExpression: number;
}

export interface Gene {
  id: string;
  name: string;
  style?: CSSProperties;
}

export interface Filters {
  datasets?: DefaultMenuSelectOption[];
  developmentStages?: DefaultMenuSelectOption[];
  diseases?: DefaultMenuSelectOption[];
  ethnicities?: DefaultMenuSelectOption[];
  sexes?: DefaultMenuSelectOption[];
}
