import { DefaultMenuSelectOption } from "@czi-sds/components";
import { CSSProperties } from "react";
import { ChartFormat } from "../components/HeatMap/components/Chart/types";
import { CellTypeMetadata } from "../../WheresMyGene/components/HeatMap/utils";

export interface Organism {
  id: string;
  name: string;
}

/** tissue name */
export type Tissue = string;

/** This object holds a cell type and its gene expressions. */
export interface CellTypeSummary {
  id: CellType["id"];
  viewId: CellType["viewId"];
  name: CellType["name"];
  total_count: number;
  isAggregated: boolean;
  geneExpressions?: {
    [geneName: string]: CellTypeGeneExpressionSummaryData;
  };
  order: CellType["order"];
}

export interface GeneInfo {
  ncbi_url: string;
  name: string;
  synonyms: string[];
  summary: string;
  show_warning_banner: boolean;
}

export interface CellType {
  viewId: ViewId;
  // ontology term id
  id: string;
  name: string;
  order: number;
  total_count: number;
  isAggregated: boolean;
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
  /** percentage of the current subset of cells that express this gene */
  pc: number;
  /** mean expression of the current subset of cells that express this gene */
  me: number;
  /** Expressed cell count */
  n: number;
  /** Tissue Composition */
  tpc: number;
}

export type CellTypeId = string;
export type CompareOptionId = string;

/**
 * (thuang): This ID is used to uniquely identify a cell type and its compare option
 * cellTypeId$compareOptionId
 * E.g., "CL:0000003$PATO:0000383", "CL:0000003$aggregated"
 **/
export type ViewId = `${CellTypeId}$${CompareOptionId}`;

/**
 * This object describes the cell type id and its gene expression metadata given
 * the current subset of cells.
 */
export interface CellTypeGeneExpressionSummaryData {
  viewId: ViewId;
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
  publications?: DefaultMenuSelectOption[];
  sexes?: DefaultMenuSelectOption[];
  tissues?: DefaultMenuSelectOption[];
}

export enum SORT_BY {
  CELL_ONTOLOGY = "CELL_ONTOLOGY",
  H_CLUSTER = "H_CLUSTER",
  USER_ENTERED = "USER_ENTERED",
  COLOR_SCALED = "SCALED",
  COLOR_UNSCALED = "UNSCALED",
}

export interface ChartProps {
  chartData: ChartFormat[];
  geneNames: string[];
  cellTypeMetadata: CellTypeMetadata[];
}
