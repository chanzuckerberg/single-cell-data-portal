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
  total_count: number;
  geneExpressions?: {
    [geneName: string]: CellTypeGeneExpressionSummaryData;
  };
}

export interface CellType {
  id: string;
  name: string;
  depth?: number;
  total_count: number;
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
  /** Expressed cell count */
  n: number;
  /** Tissue Composition */
  tpc: number;
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
  /** Tissue Composition */
  tissuePercentage: number;
  /** Expressed cell count */
  expressedCellCount: number;
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

export enum SORT_BY {
  CELL_ONTOLOGY = "CELL_ONTOLOGY",
  H_CLUSTER = "H_CLUSTER",
  USER_ENTERED = "USER_ENTERED",
  COLOR_SCALED = "SCALED",
  COLOR_UNSCALED = "UNSCALED",
}
