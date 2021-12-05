import { ChangeEvent, Dispatch, SetStateAction } from "react";
import { CellValue, Row } from "react-table";
import { Collection, IS_PRIMARY_DATA, Ontology } from "src/common/entities";
import { CategoryValueView } from "src/common/hooks/useCategoryFilter";

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

/* Filterable values of datasets and collections. */
export interface Categories {
  assay: Ontology[];
  cell_type: Ontology[];
  disease: Ontology[];
  is_primary_data: IS_PRIMARY_DATA[];
  organism: Ontology[];
  sex: Ontology[];
  tissue: Ontology[];
}

/* Join of collection, dataset and aggregated dataset information, optimized for filtering collections (that is,
   datasets grouped by collection. */
export interface CollectionRow extends Categories {
  id: string;
  name: string;
  published_at: number;
  revised_at?: number;
}

/* Join of dataset and collection information, optimized for filtering datasets */
export interface DatasetRow extends Categories {
  cell_count: number | null;
  collection_id: Collection["id"];
  collection_name: Collection["name"];
  id: string;
  name: string;
  published_at: number;
  revised_at?: number;
}

export type FilterCategoryValuesFn = (
  values: CategoryValueView[],
  searchValue: string
) => CategoryValueView[];

export enum IS_PRIMARY_DATA_LABEL {
  "PRIMARY" = "primary",
  "SECONDARY" = "composed",
}

/* Value for displaying pluralized metadata labels, for example, "tissues" or "diseases". */
export enum PLURALIZED_METADATA_LABEL {
  ASSAY = "assays",
  CELL_TYPE = "cell types",
  DISEASE = "diseases",
  ORGANISM = "organisms",
  TISSUE = "tissues",
}

/* Display value of category labels. */
export enum CATEGORY_LABEL {
  assay = "Assay",
  cell_type = "Cell Type",
  disease = "Disease",
  is_primary_data = "Data Source",
  organism = "Organism",
  tissue = "Tissue",
  sex = "Sex",
}

/* "row" prop passed to react-table's Cell function */
export type RowPropsValue = { row: Row<CollectionRow> };

// TODO(cc)
export type SetSearchValueFn = Dispatch<SetStateAction<string>>;

// TODO(cc)
export type OnUpdateSearchValueFn = (
  changeEvent: ChangeEvent<HTMLInputElement>,
  setSearchValue: SetSearchValueFn
) => void;
