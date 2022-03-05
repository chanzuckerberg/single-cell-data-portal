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
  "BETWEEN" = "between",
  "INCLUDES_SOME" = "includesSome",
}

/**
 * Filterable metadata keys.
 */
export enum CATEGORY_KEY {
  "ASSAY" = "assay",
  "CELL_COUNT" = "cell_count",
  "CELL_TYPE" = "cell_type",
  "DEVELOPMENT_STAGE_ANCESTORS" = "development_stage_ancestors",
  "DISEASE" = "disease",
  "ETHNICITY" = "ethnicity",
  "IS_PRIMARY_DATA" = "is_primary_data",
  "MEAN_GENES_PER_CELL" = "mean_genes_per_cell",
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
    categoryKey: CATEGORY_KEY.CELL_COUNT,
    categoryType: CATEGORY_FILTER_TYPE.BETWEEN,
    multiselect: false,
  },
  {
    categoryKey: CATEGORY_KEY.CELL_TYPE,
    categoryType: CATEGORY_FILTER_TYPE.INCLUDES_SOME,
    multiselect: true,
  },
  {
    categoryKey: CATEGORY_KEY.DEVELOPMENT_STAGE_ANCESTORS,
    categoryType: CATEGORY_FILTER_TYPE.INCLUDES_SOME,
    multiselect: true,
  },
  {
    categoryKey: CATEGORY_KEY.ETHNICITY,
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
    categoryKey: CATEGORY_KEY.MEAN_GENES_PER_CELL,
    categoryType: CATEGORY_FILTER_TYPE.BETWEEN,
    multiselect: false,
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
  ethnicity: Ontology[];
  is_primary_data: IS_PRIMARY_DATA[];
  organism: Ontology[];
  sex: Ontology[];
  tissue: Ontology[];
}

/**
 * View model of category.
 */
export type CategoryView =
  | OntologyCategoryView
  | RangeCategoryView
  | SelectCategoryView;

/**
 * Join of collection, dataset and aggregated dataset information, optimized for filtering collections (that is,
 * datasets grouped by collection.
 */
export interface CollectionRow extends Categories, PublisherMetadataCategories {
  id: string;
  name: string;
  published_at: number;
  publisher_metadata?: PublisherMetadata; // Collection publication metadata
  recency: number; // Used by sort
  revised_at?: number;
  summaryCitation: string;
}

/**
 * Category values to be used as keys. For example, "Homo sapiens" or "10X 3' v2 sequencing".
 */
export type CategoryValueKey = string;

/**
 * Join of dataset and collection information, optimized for filtering datasets.
 */
export interface DatasetRow extends Categories, PublisherMetadataCategories {
  cell_count: number | null;
  collection_id: Collection["id"];
  collection_name: Collection["name"];
  development_stage_ancestors: string[];
  explorer_url: string;
  id: string;
  isOverMaxCellCount: boolean;
  mean_genes_per_cell: number | null;
  name: string;
  published_at: number;
  recency: number; // Used by sort
  revised_at?: number;
}

/**
 * Display values of unspecified ethnicity labels.
 */
export enum ETHNICITY_UNSPECIFIED_LABEL {
  "na" = "N/A",
  "unknown" = "Unknown",
}

/**
 * Function returns filtered category values when category key contains filter category input value.
 */
export type FilterCategoryValuesFn = (
  values: SelectCategoryValueView[],
  searchValue: string
) => SelectCategoryValueView[];

/**
 * Function returns filtered category values with a count greater than zero.
 */
export type FilterCategoryValuesWithCountFn = (
  values: SelectCategoryValueView[]
) => SelectCategoryValueView[];

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
  | CATEGORY_KEY.CELL_COUNT
  | CATEGORY_KEY.DEVELOPMENT_STAGE_ANCESTORS
  | CATEGORY_KEY.MEAN_GENES_PER_CELL
  | CATEGORY_KEY.IS_PRIMARY_DATA
  | CATEGORY_KEY.PUBLICATION_DATE_VALUES
  | CATEGORY_KEY.PUBLICATION_AUTHORS
>;

/**
 * Node of development stage tree including ontology ID, label and any possible children.
 */
export interface DevelopmentStageNode extends Ontology {
  children?: DevelopmentStageNode[];
}

/**
 * Development stage tree structure, keyed by organism.
 */
export interface DevelopmentStageTree {
  [key: string]: DevelopmentStageNode[];
}

/**
 * Homo sapiens, Mus musculus and other organisms development stage ontology tree.
 */
/* eslint-disable sort-keys -- disabling key order for readability. */
export const DEVELOPMENT_STAGES: DevelopmentStageTree = {
  // TODO(cc) rename
  "Homo sapiens": [
    {
      label: "Prenatal (conception–birth)",
      ontology_term_id: "HsapDv:0000045",
      children: [
        {
          label: "Embryonic human (0–56 days)",
          ontology_term_id: "HsapDv:0000002",
          children: [
            {
              label: "Carnegie (CS1)",
              ontology_term_id: "HsapDv:0000003",
            },
            {
              label: "Cleavage (CS2)",
              ontology_term_id: "HsapDv:0000004",
            },
            {
              label: "Blastula (CS3–5)",
              ontology_term_id: "HsapDv:0000006",
            },
            {
              label: "Gastrula (CS6)",
              ontology_term_id: "HsapDv:0000010",
            },
            {
              label: "Neurula (CS7–8)",
              ontology_term_id: "HsapDv:0000012",
            },
            {
              label: "Organogenesis (CS9–23)",
              ontology_term_id: "HsapDv:0000015",
            },
          ],
        },
        {
          label: "Fetal (>56 days–birth)",
          ontology_term_id: "HsapDv:0000037",
        },
      ],
    },
    {
      label: "Immature (0–12 years)",
      ontology_term_id: "HsapDv:0000080",
      children: [
        {
          label: "Newborn human (0–1 month)",
          ontology_term_id: "HsapDv:0000082",
        },
        {
          label: "Infant (1–23 months)",
          ontology_term_id: "HsapDv:0000083",
        },
        {
          label: "Child (2–12 years)",
          ontology_term_id: "HsapDv:0000081",
        },
      ],
    },
    {
      label: "Mature (13+ years)",
      ontology_term_id: "HsapDv:0000204",
      children: [
        {
          label: "Adolescent (13–19 years)",
          ontology_term_id: "HsapDv:0000086",
        },
        {
          label: "Human adult (19+ years)",
          ontology_term_id: "HsapDv:0000087",
          children: [
            {
              label: "Early adulthood (19–45 years)",
              ontology_term_id: "HsapDv:0000088",
            },
            {
              label: "Late adulthood (45+ years)",
              ontology_term_id: "HsapDv:0000091",
            },
          ],
        },
      ],
    },
  ],
  "Mus musculus": [
    {
      label: "Prenatal",
      ontology_term_id: "MmusDv:0000042",
      children: [
        {
          label: "Embryonic mouse",
          ontology_term_id: "MmusDv:0000002",
          children: [
            {
              label: "Thelier stage 1 (TS1)",
              ontology_term_id: "MmusDv:0000003",
            },
            {
              label: "Cleavage (TS2–3)",
              ontology_term_id: "MmusDv:0000004",
            },
            {
              label: "Blastula (TS4–8)",
              ontology_term_id: "MmusDv:0000007",
            },
            {
              label: "Gastrula (TS9–10)",
              ontology_term_id: "MmusDv:0000013",
            },
            {
              label: "Thelier stage 11 (TS11)",
              ontology_term_id: "MmusDv:0000017",
            },
            {
              label: "Organogenesis (TS11–22)",
              ontology_term_id: "MmusDv:0000018",
            },
          ],
        },
        {
          label: "Fetal (TS23–26)",
          ontology_term_id: "MmusDv:0000031",
        },
      ],
    },
    {
      label: "Post-partum (Birth+)",
      ontology_term_id: "MmusDv:0000092",
      children: [
        {
          label: "Immature (0–6 weeks)",
          ontology_term_id: "MmusDv:0000043",
          children: [
            {
              label: "Thelier stage 27 (0–3 days)",
              ontology_term_id: "MmusDv:0000036",
            },
            {
              label: "Premature (3 days–6 weeks)",
              ontology_term_id: "MmusDv:0000112",
            },
          ],
        },
        {
          label: "Mature (6+ weeks)",
          ontology_term_id: "MmusDv:0000110",
          children: [
            {
              label: "Early adulthood (6 weeks–7 months)",
              ontology_term_id: "MmusDv:0000061",
            },
            {
              label: "Late adulthood (7+ months)",
              ontology_term_id: "MmusDv:0000097",
            },
          ],
        },
      ],
    },
  ],
  "Other organisms": [
    {
      label: "Embryo",
      ontology_term_id: "UBERON:0000068",
      children: [
        {
          label: "Zygote",
          ontology_term_id: "UBERON:0000106",
        },
        {
          label: "Cleavage",
          ontology_term_id: "UBERON:0000107",
        },
        {
          label: "Blastula",
          ontology_term_id: "UBERON:0000108",
        },
        {
          label: "Gastrula",
          ontology_term_id: "UBERON:0000109",
        },
        {
          label: "Neurula",
          ontology_term_id: "UBERON:0000110",
        },
        {
          label: "Organogenesis",
          ontology_term_id: "UBERON:0000111",
        },
        {
          label: "Late embryonic",
          ontology_term_id: "UBERON:0007220",
        },
      ],
    },
    {
      label: "Post embryonic",
      ontology_term_id: "UBERON:0000092",
      children: [
        {
          label: "Larval",
          ontology_term_id: "UBERON:0000069",
        },
        {
          label: "Pupal",
          ontology_term_id: "UBERON:0000070",
        },
        {
          label: "Nursing",
          ontology_term_id: "UBERON:0018685",
        },
        {
          label: "Fully formed",
          ontology_term_id: "UBERON:0000066",
          children: [
            {
              label: "Sexually immature",
              ontology_term_id: "UBERON:0000112",
            },
            {
              label: "Post-juvenile adult",
              ontology_term_id: "UBERON:0000113",
            },
          ],
        },
      ],
    },
  ],
};
/* eslint-enable sort-keys -- disabling key order for readability. */

/**
 * Display value of category labels.
 */
export enum CATEGORY_LABEL {
  assay = "Assay",
  cell_count = "Cell Count",
  cell_type = "Cell Type",
  development_stage_ancestors = "Development Stage",
  disease = "Disease",
  ethnicity = "Ethnicity",
  is_primary_data = "Data Source",
  mean_genes_per_cell = "Gene Count",
  publicationAuthors = "Authors",
  publicationDateValues = "Publication Date",
  organism = "Organism",
  tissue = "Tissue",
  sex = "Sex",
}

/**
 * Function invoked when selected state of a category value is toggled or range is selected.
 */
export type OnFilterFn = (
  categoryKey: CategoryKey,
  selectedValue: CategoryValueKey | Range
) => void;

/**
 * Function invoked when filter category input value is changed.
 */
export type OnUpdateSearchValueFn = (
  changeEvent: ChangeEvent<HTMLInputElement>,
  setSearchValue: SetSearchValueFn
) => void;

/**
 * Min and max values selected in range category. Empty array if no range is specified (e.g. on clear of range).
 */
export type Range = [number, number] | [];

/**
 * View model of species-specific ontology tree.
 */
export interface OntologyCategorySpeciesView {
  children: OntologyCategoryValueView[]; // TODO(cc) rename to development type?
  label: string;
}

/**
 * View model of category value, selected state, children, children selected state and count for ontology categories.
 */
export interface OntologyCategoryValueView extends SelectCategoryValueView {
  children?: OntologyCategoryValueView[];
  selectedTag?: boolean; // True if value, if selected, is to be included in the selected set of tags.
  partialSelected: boolean; // True if value is a parent node and some children node are selected. // TODO(cc) selectedPartial?
}

/**
 * View model of ontology category.
 */
export interface OntologyCategoryView {
  key: CATEGORY_KEY;
  label: CATEGORY_LABEL;
  species: OntologyCategorySpeciesView[];
}

/**
 * View model of range metadata key.
 */
export interface RangeCategoryView {
  key: CATEGORY_KEY;
  label: CATEGORY_LABEL;
  max: number;
  min: number;
  selectedMax?: number;
  selectedMin?: number;
}

/**
 * Possible set of values that publication dates can be binned into.
 */
export const PUBLICATION_DATE_VALUES: number[] = [1, 3, 6, 12, 24, 36];

/**
 * Display values of publicationDateValue labels. Enum must be non-numeric key for reverse lookup.
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
 * View model of metadata value, selected state and count for single or multiselect categories.
 */
export interface SelectCategoryValueView {
  key: CategoryValueKey;
  count: number;
  label: string;
  selected: boolean;
}

/**
 * Metadata values grouped by metadata key, for single or multiselect categories.
 */
export interface SelectCategoryView {
  key: CATEGORY_KEY;
  label: CATEGORY_LABEL;
  values: SelectCategoryValueView[];
}

/**
 * Function invoked to update state for the filter category input value.
 */
export type SetSearchValueFn = Dispatch<SetStateAction<string>>;
