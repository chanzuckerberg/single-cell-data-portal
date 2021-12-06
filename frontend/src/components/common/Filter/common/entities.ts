import { ChangeEvent, Dispatch, SetStateAction } from "react";
import { CellValue, Row } from "react-table";
import { Collection, IS_PRIMARY_DATA, Ontology } from "src/common/entities";
import { CategoryKey } from "src/common/hooks/useCategoryFilter";

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

/* Metadata values grouped by metadata key. */
export interface CategoryView {
  key: CATEGORY_KEY;
  label: CATEGORY_LABEL;
  values: CategoryValueView[];
}

/* "value" prop passed to react-table's Cell function */
export type CellPropsValue<T> = { value: CellValue<T> };

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

/* Category values to be used as keys. For example, "Homo sapiens" or "10X 3' v2 sequencing". */
export type CategoryValueKey = string;

/* View model of metadata value, selected state and count. */
export interface CategoryValueView {
  count: number;
  key: CategoryValueKey;
  label: string;
  selected: boolean;
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

/* Filterable metadata keys where the type of the corresponding value is Ontology. Currently, that is all metadata
   keys except is_primary_data. */
export type OntologyCategoryKey = keyof Omit<
  Record<CATEGORY_KEY, string>,
  CATEGORY_KEY.IS_PRIMARY_DATA
>;

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

/* Function invoked when selected state of a category value is toggled. */
export type OnFilterFn = (
  categoryKey: CategoryKey,
  categoryValueKey: CategoryValueKey
) => void;

/* Function invoked when filter category input value is changed. */
export type OnUpdateSearchValueFn = (
  changeEvent: ChangeEvent<HTMLInputElement>,
  setSearchValue: SetSearchValueFn
) => void;

/* "row" prop passed to react-table's Cell function. */
export type RowPropsValue<T extends Categories> = { row: Row<T> };

/* Function invoked to update state for the filter category input value. */
export type SetSearchValueFn = Dispatch<SetStateAction<string>>;
