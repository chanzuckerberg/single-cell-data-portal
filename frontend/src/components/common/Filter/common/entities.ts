import { ChangeEvent, Dispatch, SetStateAction } from "react";
import { CellValue, Row } from "react-table";
import {
  Collection,
  IS_PRIMARY_DATA,
  Ontology,
  PublisherMetadata,
} from "src/common/entities";
import { CategoryKey } from "src/common/hooks/useCategoryFilter";

/**
 * Configuration model of category.
 */
export interface CategoryConfig {
  categoryKey: CATEGORY_KEY;
  categoryType: CATEGORY_FILTER_TYPE;
  multiselect: boolean; // True if category can have multiple values selected.
}

/**
 * Possible types of category filters, matches React Table's filter types.
 */
export enum CATEGORY_FILTER_TYPE {
  "INCLUDES_SOME" = "includesSome",
}

/**
 * Filterable metadata keys.
 */
export enum CATEGORY_KEY {
  "ASSAY" = "assay",
  "CELL_TYPE" = "cell_type",
  "DISEASE" = "disease",
  "IS_PRIMARY_DATA" = "is_primary_data",
  "ORGANISM" = "organism",
  "PUBLICATION_AUTHORS" = "publicationAuthors",
  "PUBLICATION_DATE_VALUES" = "publicationDateValues",
  "SEX" = "sex",
  "TISSUE" = "tissue",
}

/**
 * Configuration for each category.
 */
const CATEGORY_CONFIGS: CategoryConfig[] = [
  {
    categoryKey: CATEGORY_KEY.ASSAY,
    categoryType: CATEGORY_FILTER_TYPE.INCLUDES_SOME,
    multiselect: true,
  },
  {
    categoryKey: CATEGORY_KEY.CELL_TYPE,
    categoryType: CATEGORY_FILTER_TYPE.INCLUDES_SOME,
    multiselect: true,
  },
  {
    categoryKey: CATEGORY_KEY.DISEASE,
    categoryType: CATEGORY_FILTER_TYPE.INCLUDES_SOME,
    multiselect: true,
  },
  {
    categoryKey: CATEGORY_KEY.IS_PRIMARY_DATA,
    categoryType: CATEGORY_FILTER_TYPE.INCLUDES_SOME,
    multiselect: true,
  },
  {
    categoryKey: CATEGORY_KEY.ORGANISM,
    categoryType: CATEGORY_FILTER_TYPE.INCLUDES_SOME,
    multiselect: true,
  },
  {
    categoryKey: CATEGORY_KEY.PUBLICATION_AUTHORS,
    categoryType: CATEGORY_FILTER_TYPE.INCLUDES_SOME,
    multiselect: true,
  },
  {
    categoryKey: CATEGORY_KEY.PUBLICATION_DATE_VALUES,
    categoryType: CATEGORY_FILTER_TYPE.INCLUDES_SOME,
    multiselect: false,
  },
  {
    categoryKey: CATEGORY_KEY.SEX,
    categoryType: CATEGORY_FILTER_TYPE.INCLUDES_SOME,
    multiselect: true,
  },
  {
    categoryKey: CATEGORY_KEY.TISSUE,
    categoryType: CATEGORY_FILTER_TYPE.INCLUDES_SOME,
    multiselect: true,
  },
];

/**
 * Category configs keyed by category key, for convenience. Using object literal with type KeyedCategoryConfig rather
 * than generic Map to prevent having to null check values.
 */
export const CATEGORY_CONFIGS_BY_CATEGORY_KEY: KeyedCategoryConfigs =
  CATEGORY_CONFIGS.reduce(
    (accum: KeyedCategoryConfigs, config: CategoryConfig) => {
      return {
        ...accum,
        [config.categoryKey]: config,
      };
    },
    {} as KeyedCategoryConfigs
  );

/**
 * Metadata values grouped by metadata key.
 */
export interface CategoryView {
  key: CATEGORY_KEY;
  label: CATEGORY_LABEL;
  values: CategoryValueView[];
}

/**
 * "value" prop passed to react-table's Cell function
 */
export type CellPropsValue<T> = { value: CellValue<T> };

/**
 * Filterable values of datasets and collections.
 */
export interface Categories {
  assay: Ontology[];
  cell_type: Ontology[];
  disease: Ontology[];
  is_primary_data: IS_PRIMARY_DATA[];
  organism: Ontology[];
  sex: Ontology[];
  tissue: Ontology[];
}

/**
 * Join of collection, dataset and aggregated dataset information, optimized for filtering collections (that is,
 * datasets grouped by collection.
 */
export interface CollectionRow extends Categories, PublisherMetadataCategories {
  id: string;
  name: string;
  published_at: number;
  publisher_metadata?: PublisherMetadata; // Collection publication metadata
  revised_at?: number;
}

/**
 * Category values to be used as keys. For example, "Homo sapiens" or "10X 3' v2 sequencing".
 */
export type CategoryValueKey = string;

/**
 * View model of metadata value, selected state and count.
 */
export interface CategoryValueView {
  count: number;
  key: CategoryValueKey;
  label: string;
  selected: boolean;
}

/**
 * Join of dataset and collection information, optimized for filtering datasets.
 */
export interface DatasetRow extends Categories, PublisherMetadataCategories {
  cell_count: number | null;
  collection_id: Collection["id"];
  collection_name: Collection["name"];
  explorer_url: string;
  id: string;
  isOverMaxCellCount: boolean;
  name: string;
  publisher_metadata?: PublisherMetadata; // TODO(cc) remove before PR, required for temp col display values
  published_at: number;
  revised_at?: number;
}

/**
 * Function returns filtered category values when category key contains filter category input value.
 */
export type FilterCategoryValuesFn = (
  values: CategoryValueView[],
  searchValue: string
) => CategoryValueView[];

/**
 * Function returns filtered category values with a count greater than zero.
 */
export type FilterCategoryValuesWithCountFn = (
  values: CategoryValueView[]
) => CategoryValueView[];

/**
 * Model of category configs keyed by category key. Used instead of generic Map to prevent null checking when grabbing
 * keyed value.
 */
type KeyedCategoryConfigs = { [K in CATEGORY_KEY]: CategoryConfig };

/**
 * Display values of is_primary_data labels.
 */
export enum IS_PRIMARY_DATA_LABEL {
  "PRIMARY" = "primary",
  "SECONDARY" = "composed",
}

/**
 * Filterable metadata keys where the type of the corresponding value is Ontology. Currently, that is all metadata
 * keys except is_primary_data.
 */
export type OntologyCategoryKey = keyof Omit<
  Record<CATEGORY_KEY, string>,
  | CATEGORY_KEY.IS_PRIMARY_DATA
  | CATEGORY_KEY.PUBLICATION_DATE_VALUES
  | CATEGORY_KEY.PUBLICATION_AUTHORS
>;

/**
 * Display value of category labels.
 */
export enum CATEGORY_LABEL {
  assay = "Assay",
  cell_type = "Cell Type",
  disease = "Disease",
  is_primary_data = "Data Source",
  publicationAuthors = "Authors",
  publicationDateValues = "Publication Date",
  organism = "Organism",
  tissue = "Tissue",
  sex = "Sex",
}

/**
 * Function invoked when selected state of a category value is toggled.
 */
export type OnFilterFn = (
  categoryKey: CategoryKey,
  categoryValueKey: CategoryValueKey
) => void;

/**
 * Function invoked when filter category input value is changed.
 */
export type OnUpdateSearchValueFn = (
  changeEvent: ChangeEvent<HTMLInputElement>,
  setSearchValue: SetSearchValueFn
) => void;

/**
 * Possible set of values that publication dates can be binned into.
 */
export const PUBLICATION_DATE_VALUES: number[] = [1, 3, 6, 12, 24, 36];

/**
 * Display values of publicationDateValue labels. Note, enum must be non-numeric key for reverse lookup.
 */
export enum PUBLICATION_DATE_LABELS {
  "LABEL_1" = "1 month",
  "LABEL_12" = "12 months",
  "LABEL_24" = "2 years",
  "LABEL_3" = "3 months",
  "LABEL_36" = "3 years",
  "LABEL_6" = "6 months",
}

/**
 * Publication-related filterable values of collections and datasets.
 */
export interface PublisherMetadataCategories {
  publicationAuthors?: string[];
  publicationDateValues?: number[]; // Set of date bins that publication date falls within
}

/**
 * "row" prop passed to react-table's Cell function.
 */
export type RowPropsValue<T extends Categories> = { row: Row<T> };

/**
 * Function invoked to update state for the filter category input value.
 */
export type SetSearchValueFn = Dispatch<SetStateAction<string>>;
