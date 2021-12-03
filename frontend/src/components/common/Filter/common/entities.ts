import { CellValue, Row } from "react-table";
import { Collection, IS_PRIMARY_DATA, Ontology } from "src/common/entities";

/* Filterable metadata keys */
export enum CATEGORY_KEY {
  "ASSAY" = "assay",
  "CELL_TYPE" = "cell_type",
  "DISEASE" = "disease",
  "IS_PRIMARY_DATA" = "is_primary_data",
  "ORGANISM" = "organism",
  "SEX" = "sex",
  "TISSUE" = "tissue",
}

/* Filterable metadata object key. For example, "assay", "cell_type" or "is_primary_data". Used for object key lookups */
export type CategoryKey = keyof Record<CATEGORY_KEY, string>;

/* "value" prop passed to react-table's Cell function */
export type CellPropsValue = { value: CellValue<string[]> };

export interface CollectionCategoryValues {
  assayAggregated: Ontology[];
  cellTypeAggregated: Ontology[];
  diseaseAggregated: Ontology[];
  isPrimaryDataAggregated: IS_PRIMARY_DATA[];
  organismAggregated: Ontology[];
  sexAggregated: Ontology[];
  tissueAggregated: Ontology[];
}

/* Join of collection, dataset and aggregated dataset information, optimized for filtering collections (that is,
   datasets grouped by collection. */
export type CollectionRow = DatasetRow & CollectionCategoryValues;

/* Join of dataset and collection information, optimized for filtering datasets */
export interface DatasetRow {
  assay: Ontology[];
  cell_count: number | null;
  cell_type: Ontology[];
  collection_id: Collection["id"];
  collection_name: Collection["name"];
  disease: Ontology[];
  id: string;
  is_primary_data: IS_PRIMARY_DATA[]; // Handle "BOTH" as ["primary", "secondary"]
  name: string;
  organism: Ontology[];
  published_at: number;
  revised_at?: number;
  sex: Ontology[];
  tissue: Ontology[];
}

/* "row" prop passed to react-table's Cell function */
export type RowPropsValue = { row: Row<CollectionRow> };
