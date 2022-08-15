import { ChangeEvent, Dispatch, SetStateAction } from "react";
import { CellValue, Row } from "react-table";
import { EVENTS } from "src/common/analytics/events";
import { Collection, Ontology, PublisherMetadata } from "src/common/entities";
import { CategoryFilterId } from "src/common/hooks/useCategoryFilter";

// ** Label discriminating unions ** //

/**
 * Model of a filter that uses an ID value for filtering and must perform a lookup for display.
 */
export interface LookupByTermIdLabelKind {
  labelKind: "LOOKUP_LABEL_BY_TERM_ID";
}

/**
 * Model of a filter that uses the same value for both filtering and display. For example, publication date, cell count,
 * assay etc.
 */
export interface ValueLabelKind {
  labelKind: "VALUE";
}

// ** Value source kind discriminating unions ** //

/**
 * Model of a filter that uses a hand-curated set of terms to determine the filter category values allowed for display.
 * For example, development stage, tissue system and tissue organ.
 */
export interface CuratedValueSourceKind {
  // TODO(cc) - revisit booleans below, can we remove?
  isLabelVisible: boolean; // True if ontology label is to be displayed (e.g. "Other Organisms")
  isSearchable: boolean; // True if ontology category search box is displayed
  isZerosVisible: boolean; // True if zero values are to be displayed
  mask: OntologyTermSet; // TODO(cc) rename mask to valueSource? rename OntologyTermSet
  valueSourceKind: "CURATED";
}

/**
 * Model of a filter that has no restrictions on the filter category values allowed for display. For example, assay,
 * publication date etc.
 */
export interface NoneValueSourceKind {
  valueSourceKind: "NONE";
}

// ** Match discriminating unions ** //

/**
 * Model of a filter that uses a "between" method for filtering values. Terminology matches React Table filters. Used
 * by filter categories that are ranges.
 */
export interface BetweenMatchKind {
  matchKind: "BETWEEN";
}

/**
 * Model of a filter that uses an "includes some" method for filtering values. Terminology matches React Table filters.
 * Used by any filter category that has checkboxes, either straight-up select filters or ontology-backed filters.
 */
export interface IncludesSomeMatchKind {
  matchKind: "INCLUDES_SOME";
}

// ** Query discriminating unions ** //

/**
 * Model of a filter that uses a standard faceted search paradigm when determining the set of queries to apply to itself
 * on filter. That is, it uses all currently selected values - except its own selected values - to filter. For example,
 * publication date, assay etc.
 */
export interface ExcludesSelfQueryKind {
  queryKind: "EXCLUDES_SELF";
}

/**
 * Model of a filter that uses a modified faceted search paradigm when determining the set of queries to apply to
 * itself on filter: it ignores its own selected values as well as selected values in its children category filters. For
 * example, the tissue system category filter should ignore selected tissue systems values as well as selected values
 * in tissue organ and tissue.
 */
export interface ExcludesSelfAndChildrenQueryKind {
  childrenCategoryFilterIds: CATEGORY_FILTER_ID[];
  queryKind: "EXCLUDES_SELF_AND_CHILDREN";
}

// ** Value restriction discriminating unions ** //

/**
 * Model of a filter that uses the set of selected values in a parent filter to determine the filter category values
 * allowed for display. For example, tissue. NOTE! The order of parents defined in parentCategoryFilterIds drives
 * the order of the restrictions applied to the category filter. For example, for tissue, organ must be specified
 * before system so that when restricting tissue values, if both system and organ have a selected value, then only the
 * selected organ is applied as a restriction. If both selected system and selected organ were applied to tissue, the
 * tissue values would over-show.
 */
export interface SelectedParentTermsValueRestrictionKind {
  valueRestrictionKind: "CHILDREN_OF_SELECTED_PARENT_TERMS";
  parentCategoryFilterIds: CATEGORY_FILTER_ID[];
}

/**
 * Model of a filter that has no restrictions on the filter category values allowed for display. For example, assay,
 * publication date etc.
 */
export interface NoneValueRestrictionKind {
  valueRestrictionKind: "NONE";
}

// ** Category filter constants ** //

/**
 * Complete set of category filter IDs that can possibly be included in filter.
 */
export enum CATEGORY_FILTER_ID { // TODO(cc) rename usage (ie variables to match new name here)
  "ASSAY" = "ASSAY",
  "CELL_COUNT" = "CELL_COUNT",
  "CELL_TYPE_DEPRECATED" = "CELL_TYPE_DEPRECATED",
  "CELL_TYPE" = "CELL_TYPE",
  "DEVELOPMENT_STAGE" = "DEVELOPMENT_STAGE",
  "DISEASE" = "DISEASE",
  "ETHNICITY" = "ETHNICITY",
  "GENE_COUNT" = "GENE_COUNT",
  "ORGANISM" = "ORGANISM",
  "PUBLICATION_AUTHORS" = "PUBLICATION_AUTHORS",
  "PUBLICATION_DATE_VALUES" = "PUBLICATION_DATE_VALUES",
  "SEX" = "SEX",
  "TISSUE_DEPRECATED" = "TISSUE_DEPRECATED", // TODO(cc) remove with #2569.
  "TISSUE" = "TISSUE",
  "TISSUE_ORGAN" = "TISSUE_ORGAN",
  "TISSUE_SYSTEM" = "TISSUE_SYSTEM",
}

// ** Category filter types ** //

/**
 * Base category filter configuration model containing values shared across the different types of category filters.
 */
export interface BaseCategoryFilterConfig {
  analyticsEvent?: EVENTS;
  categoryFilterId: CATEGORY_FILTER_ID;
  filterOnKey: FilterKey; // Key in result set row values to filter on.
  label: string;
  multiselect: boolean; // True if category can have multiple values selected.
  pinnedCategoryValues?: CATEGORY_VALUE_KEY[];
  tooltip?: string;
}

/**
 * Filter category that uses checkboxes for filtering, is ontology-aware, and uses a hand-curated set of ontology terms
 * (or "mask") to determine the values to display as filter category values. For example, development stage, tissue
 * system and tissue organ.
 */
export type CuratedOntologyCategoryFilterConfig = BaseCategoryFilterConfig &
  IncludesSomeMatchKind &
  LookupByTermIdLabelKind &
  CuratedValueSourceKind &
  ExcludesSelfQueryKind &
  NoneValueRestrictionKind;

/**
 * Filter category that uses checkboxes for filtering, is ontology-aware, and whose displayable values are determined
 * by parent category filters. For example, tissue.
 */
export type LeafOntologyCategoryFilterConfig = BaseCategoryFilterConfig &
  IncludesSomeMatchKind &
  LookupByTermIdLabelKind &
  NoneValueSourceKind &
  ExcludesSelfQueryKind &
  SelectedParentTermsValueRestrictionKind;

/**
 * Filter category that uses checkboxes for filtering, is ontology-aware, and whose selected values restrict allowed
 * displayable values in children category filters. For example, tissue system.
 */
export type ParentOntologyCategoryFilterConfig = BaseCategoryFilterConfig &
  IncludesSomeMatchKind &
  LookupByTermIdLabelKind &
  CuratedValueSourceKind &
  ExcludesSelfAndChildrenQueryKind &
  NoneValueRestrictionKind;

/**
 * Filter category that uses checkboxes for filtering, is ontology-aware, whose selected values restrict allowed
 * displayable values in children category filters and whose displayable values are determined by parent category
 * filters. For example, tissue organ.
 */
export type ChildOntologyCategoryFilterConfig = BaseCategoryFilterConfig &
  IncludesSomeMatchKind &
  LookupByTermIdLabelKind &
  CuratedValueSourceKind &
  ExcludesSelfAndChildrenQueryKind &
  SelectedParentTermsValueRestrictionKind;

/**
 * Filter category that a range slider for filtering. For example, cell count and gene count.
 */
export type RangeCategoryFilterConfig = BaseCategoryFilterConfig &
  BetweenMatchKind &
  ValueLabelKind &
  NoneValueSourceKind &
  ExcludesSelfQueryKind &
  NoneValueRestrictionKind;

/**
 * Filter category that uses checkboxes for filtering and is not ontology-aware. For example, publication date and
 * author. Currently, ontology values such as assay, organism, disease etc. are also this type of filter category as,
 * although they are based on an ontology term, they are not ontology-aware.
 */
export type SelectCategoryFilterConfig = BaseCategoryFilterConfig &
  IncludesSomeMatchKind &
  ValueLabelKind &
  NoneValueSourceKind &
  ExcludesSelfQueryKind &
  NoneValueRestrictionKind;

/**
 * Set of all possible filter category configuration types.
 */
export type CategoryFilterConfig =
  | CuratedOntologyCategoryFilterConfig
  | LeafOntologyCategoryFilterConfig
  | ParentOntologyCategoryFilterConfig
  | ChildOntologyCategoryFilterConfig
  | RangeCategoryFilterConfig
  | SelectCategoryFilterConfig;

/**
 * UI configuration for each category.
 */
export interface CategoryFilterUIConfig {
  categoryFilterConfigIds: CATEGORY_FILTER_ID[];
  label: string;
}

/**
 * Possible set of keys to filter over.
 */
export type FilterKey = keyof CollectionRow | keyof DatasetRow;

/**
 * Category value keys.
 */
export enum CATEGORY_VALUE_KEY {
  NORMAL = "normal",
}

/**
 * Filterable values of datasets and collections.
 */
export interface Categories {
  assay: Ontology[];
  cell_type: Ontology[];
  cell_type_ancestors: string[];
  disease: Ontology[];
  development_stage_ancestors: string[];
  ethnicity: Ontology[];
  organism: Ontology[];
  sex: Ontology[];
  tissue: Ontology[]; // TODO(cc) remove with #2569.
  tissue_ancestors: string[];
}

/**
 * Keys of categories that are ontology arrays.
 */
export type CategoriesKeyOfTypeOntologyArray = {
  [K in keyof Categories]: Categories[K] extends Ontology[] ? K : never;
}[keyof Categories];

/**
 * View model of category view container, possibly containing more than one category to display in multiple filter
 * panels (e.g. tissue).
 */
export interface CategoryViews {
  categoryViews: CategoryView[];
  label: string;
}

/**
 * View model of category.
 */
export type CategoryView =
  | OntologyCategoryView
  | RangeCategoryView
  | SelectCategoryView;

/**
 * Category values to be used as keys. For example, "Homo sapiens" or "10X 3' v2 sequencing".
 */
export type CategoryValueKey = string;

/**
 * "value" prop passed to react-table's Cell function
 */
export type CellPropsValue<T> = { value: CellValue<T> };

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
 * Join of dataset and collection information, optimized for filtering datasets.
 */
export interface DatasetRow extends Categories, PublisherMetadataCategories {
  cell_count: number | null;
  collection_id: Collection["id"];
  collection_name: Collection["name"];
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
  "unknown" = "Unknown",
}

/**
 * List of ethnicity ontology labels to exclude from filter functionality.
 */
export const ETHNICITY_DENY_LIST = ["na"];

/**
 * Model of filter category configs keyed by category key. Used instead of generic Map to prevent null checking when
 * grabbing keyed value.
 */
export type KeyedCategoryFilterConfigs = {
  [K in CATEGORY_FILTER_ID]: CategoryFilterConfig;
};

/**
 * Possible set of organism values.
 */
export enum ORGANISM {
  "HOMO_SAPIENS" = "Homo sapiens",
  "MUS_MUSCULUS" = "Mus musculus",
}

/**
 * Function invoked when selected state of a category value is toggled or range is selected.
 */
export type OnFilterFn = (
  categoryFilterId: CategoryFilterId,
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
 * Tree view model of ontology view. For example, development stage has three tree views (human, mouse and other)
 * whereas tissue has one tree view.
 */
export interface OntologyCategoryTreeView {
  children: OntologyCategoryTreeNodeView[];
  label?: string;
  selectedViews: OntologyCategoryTreeNodeView[];
}

/**
 * View model of a node in an ontology tree view including partial selected state and possibly any children of this
 * node.
 */
export interface OntologyCategoryTreeNodeView extends SelectCategoryValueView {
  children?: OntologyCategoryTreeNodeView[];
  selectedPartial: boolean; // True if value is a parent node and some children node are selected.
}

/**
 * View model of ontology category.
 */
export interface OntologyCategoryView {
  isDisabled?: boolean;
  isSearchable: boolean;
  isZerosVisible: boolean;
  key: CategoryFilterId; // TODO(cc) rename to categoryFilterId
  label: string;
  views: OntologyCategoryTreeView[];
  tooltip?: string;
}

/**
 * Node of ontology tree including ontology ID, label and any possible children.
 */
export interface OntologyNode extends Ontology {
  children?: OntologyNode[];
}

/*
 * Ontology view keys used to lookup ontology nodes for display . For example, for the ontology ID "HsapDv:0000045",
 * the key is "HsapDv".
 */
export enum ONTOLOGY_VIEW_KEY {
  "CL" = "CL",
  "HsapDv" = "HsapDv",
  "MmusDv" = "MmusDv",
  "UBERON" = "UBERON",
}

/**
 * Labels for displaying ontology views. Currently only used by development stage filter.
 */
export enum ONTOLOGY_VIEW_LABEL {
  "HsapDv" = "Homo Sapiens",
  "MmusDv" = "Mus Musculus",
  "UBERON" = "Other Organisms",
}

/**
 * Ontology tree structures, keyed by view key. This is the allowed set of ontology values, configured per category.
 */
export type OntologyTermSet = { [K in ONTOLOGY_VIEW_KEY]?: OntologyNode[] };

/**
 * Min and max values selected in range category. Empty array if no range is specified (e.g. on clear of range).
 */
export type Range = [number, number] | [];

/**
 * View model of range metadata key.
 */
export interface RangeCategoryView {
  isDisabled?: boolean;
  key: CategoryFilterId;
  label: string;
  max: number;
  min: number;
  selectedMax?: number;
  selectedMin?: number;
  tooltip?: string;
}

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
 * Possible set of values that publication dates can be binned into.
 */
export const PUBLICATION_DATE_VALUES: number[] = [1, 3, 6, 12, 24, 36];

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
  isDisabled?: boolean;
  key: CategoryFilterId;
  label: string;
  pinnedValues: SelectCategoryValueView[];
  tooltip?: string;
  unpinnedValues: SelectCategoryValueView[];
  values: SelectCategoryValueView[]; // both pinned and unpinned values
}

/**
 * Function invoked to update state for the filter category input value.
 */
export type SetSearchValueFn = Dispatch<SetStateAction<string>>;
