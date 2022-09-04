import { CellValue, Row } from "react-table";
import { EVENTS } from "src/common/analytics/events";
import { Collection, Ontology, PublisherMetadata } from "src/common/entities";

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
 * For example, development stage.
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

// ** View kind discriminating unions ** //

/**
 * Model of filter that is ontology-aware and possibly displays multiple view panels (menus) with de/selectble
 * menu items. For example, tissue.
 */
export interface MultiPanelViewKind {
  descendants: OntologyDescendants;
  panels: CategoryFilterPanelConfig[];
  viewKind: "MULTI_PANEL";
}

/**
 * Model of filter that is ontology-aware and displays a hierarchy of ontology values possible across multiple panels.
 * For example, development stage.
 */
export interface CuratedOntologyViewKind {
  viewKind: "CURATED_ONTOLOGY";
}

/**
 * Model of filter that is not ontology-aware and displays a single view panel with a range slider. For example, cell
 * count etc.
 */
export interface RangeViewKind {
  viewKind: "RANGE";
}

/**
 * Model of filter that is not ontology-aware and displays as a single view panel (menu) with de/selectable menu items.
 * For example, assay, disease etc.
 */
export interface SelectViewKind {
  viewKind: "SELECT";
}

// ** Category filter constants ** //

/**
 * Complete set of category filter IDs that can possibly be included in filter.
 */
export enum CATEGORY_FILTER_ID {
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
  "TISSUE_CALCULATED" = "TISSUE_CALCULATED",
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
  CuratedOntologyViewKind;

/**
 * Filter category that is ontology-aware and displays a set of panels.
 * TODO(cc) rename to MultiPanelOntology  (and corresponding view too)
 */
export type OntologyMultiPanelFilterConfig = BaseCategoryFilterConfig &
  IncludesSomeMatchKind &
  LookupByTermIdLabelKind &
  NoneValueSourceKind &
  MultiPanelViewKind;

/**
 * Filter category that uses a range slider for filtering. For example, cell count and gene count.
 */
export type RangeCategoryFilterConfig = BaseCategoryFilterConfig &
  BetweenMatchKind &
  ValueLabelKind &
  NoneValueSourceKind &
  RangeViewKind;

/**
 * Filter category that uses checkboxes for filtering and is not ontology-aware. For example, publication date and
 * author. Currently, ontology values such as assay, organism, disease etc. are also this type of filter category as,
 * although they are based on an ontology term, they are not ontology-aware.
 */
export type SelectCategoryFilterConfig = BaseCategoryFilterConfig &
  IncludesSomeMatchKind &
  ValueLabelKind &
  NoneValueSourceKind &
  SelectViewKind;

/**
 * Set of all possible filter category configuration types.
 */
export type CategoryFilterConfig =
  | CuratedOntologyCategoryFilterConfig
  | OntologyMultiPanelFilterConfig
  | RangeCategoryFilterConfig
  | SelectCategoryFilterConfig;

// ** Panel value restriction discriminating unions ** //

/**
 * Model of a filter panel that uses the set of selected values in a parent filter to determine the filter panel
 * category values allowed for display. For example, the tissue panel inside the tissue filter. NOTE! The order of
 * parents defined in parentCategoryFilterIds drives the order of the restrictions applied to the filter panel. For
 * example, for tissue, organ must be specified before system so that when restricting tissue values, a selected organ
 * value can "override" a selected system value where the organ is part of the system.
 *
 */
export interface SelectedParentTermsValueRestrictionKind {
  valueRestrictionKind: "CHILDREN_OF_SELECTED_PARENT_TERMS";
  parentCategoryPanelFilterIds: CATEGORY_FILTER_PANEL_ID[];
}

/**
 * Model of a filter that has no restrictions on the filter category values allowed for display. For example, assay,
 * publication date etc.
 */
export interface NoneValueRestrictionKind {
  valueRestrictionKind: "NONE";
}

// ** Panel search discriminating unions ** //

/**
 * Model of filter that when a search term is selected, the search input is cleared. For example, system or organ.
 * TODO(cc) naming of this and below
 */
export interface SearchSingleSelectKind {
  searchKind: "SEARCH_SINGLE_SELECT";
}

/**
 * Model of filter that when a search term is selected, the search input is not cleared. For example, tissue.
 */
export interface SearchMultiSelect {
  searchKind: "SEARCH_MULTI_SELECT";
}

/**
 * Model of a filter panel that uses a hand-curated set of terms to determine the filter category values allowed for
 * display. For example, tissue system and tissue organ.
 */
interface OnlyCuratedValueSourceKind {
  mask: OntologyTermSet; // TODO(cc) rename mask to valueSource? rename OntologyTermSet same with cateogry source kind
  sourceKind: "CURATED_CATEGORIES";
}

/**
 * Model of a filter panel that displays any value not displayed in other panels. For example, tissue.
 */
interface ExceptCuratedValueSourceKind {
  sourceKind: "EXCEPT_CURATED";
}

// ** Filter value kind discriminating unions ** //

/**
 * Model of a filter panel whose category filter values filter for both inferred and explicit values. For example, tissue
 * system and tissue organ.
 */
interface InferredExplicitFilterValueKind {
  filterValueKind: "INFERRED_EXPLICIT";
}

/**
 * Model of a filter panel whose category filter values filter for only explicit values. For example, tissue.
 */
interface ExplicitOnlyFilterValueKind {
  filterValueKind: "EXPLICIT_ONLY";
}

// ** Category filter panel types ** //

/**
 * Base category filter panel configuration model containing values shared across the different types of category
 * filter panels.
 */
interface BaseCategoryFilterPanelConfig {
  id: CATEGORY_FILTER_PANEL_ID;
  label: string;
}

/**
 * Complete set of filter panel IDs.
 */
export enum CATEGORY_FILTER_PANEL_ID {
  "TISSUE_SYSTEM" = "TISSUE_SYSTEM",
  "TISSUE_ORGAN" = "TISSUE_ORGAN",
  "TISSUE" = "TISSUE",
}

/**
 * Filter category panel that uses checkboxes for filtering, is ontology-aware, and whose selected values restrict
 * allowed displayable values in children category filters. For example, tissue system.
 */
type ChildCategoryFilterPanelConfig = BaseCategoryFilterPanelConfig &
  OnlyCuratedValueSourceKind &
  InferredExplicitFilterValueKind &
  SelectedParentTermsValueRestrictionKind &
  SearchSingleSelectKind;

/**
 * Filter category panel that uses checkboxes for filtering, is ontology-aware, and whose selected values restrict
 * allowed displayable values in children category filters. For example, tissue system.
 */
type LeafCategoryFilterPanelConfig = BaseCategoryFilterPanelConfig &
  ExceptCuratedValueSourceKind &
  ExplicitOnlyFilterValueKind &
  SelectedParentTermsValueRestrictionKind &
  SearchMultiSelect;

/**
 * Filter category panel that uses checkboxes for filtering, is ontology-aware, and whose selected values restrict
 * allowed displayable values in children category filters. For example, tissue system.
 */
type ParentCategoryFilterPanelConfig = BaseCategoryFilterPanelConfig &
  OnlyCuratedValueSourceKind &
  InferredExplicitFilterValueKind &
  NoneValueRestrictionKind &
  SearchSingleSelectKind;

export type CategoryFilterPanelConfig =
  | ChildCategoryFilterPanelConfig
  | LeafCategoryFilterPanelConfig
  | ParentCategoryFilterPanelConfig;

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
  tissue: Ontology[]; // TODO(cc) revisit with #2569. don't need in catgories as we no longer filter on tissue?
  tissue_ancestors: string[]; // TODO(cc) revisit with #2569. don't need in catgories as we no longer filter on tissue_ancestors?
  tissueCalculated: string[];
}

/**
 * Keys of categories that are ontology arrays.
 */
export type CategoriesKeyOfTypeOntologyArray = {
  [K in keyof Categories]: Categories[K] extends Ontology[] ? K : never;
}[keyof Categories];

/**
 * Possible category view model types.
 */
export type CategoryView =
  | OntologyCategoryView
  | OntologyMultiPanelCategoryView
  | RangeCategoryView
  | SelectCategoryView;

/**
 * Category values to be used as keys. For example, "Homo sapiens" or "10X 3' v2 sequencing".
 * TOOD(cc) rename variables that are called categoryValueKey to categoryValueId, here and throughout
 */
export type CategoryValueId = string;

/**
 * "value" prop passed to react-table's Cell function
 */
export type CellPropsValue<T> = { value: CellValue<T> };

/**
 * Join of collection, dataset and aggregated dataset information, optimized for filtering collections (that is,
 * datasets grouped by collection.
 */
export interface CollectionRow
  extends Categories,
    PublisherMetadataCategories,
    TissueCategories {
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
export interface DatasetRow
  extends Categories,
    PublisherMetadataCategories,
    TissueCategories {
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
 * Function invoked when selected state of a category value is toggled or range is selected. Selected value is either
 * a single category value key for non ontology-aware category filters, an array of category value keys for ontology-
 * aware category filters or a range.
 */
export type OnFilterFn = (
  categoryFilterId: CATEGORY_FILTER_ID,
  categoryValueKey: CategoryValueId | null, // null for ranges. TODO(cc) remove?
  selectedValue: CategoryValueId | Range,
  source?: ON_FILTER_SOURCE
) => void;

/**
 * Location of click triggering filter action, either a value in a filter menu or the "x" in a selected tag.
 */
export enum ON_FILTER_SOURCE {
  FILTER = "FILTER",
  TAG = "TAG",
}

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
  key: CATEGORY_FILTER_ID; // TODO(cc) rename to categoryFilterId
  label: string;
  views: OntologyCategoryTreeView[];
  tooltip?: string;
}

/**
 * View model of an ontology-aware filter that contains multiple panels. For example, the overall tissue filter category.
 */
export interface OntologyMultiPanelCategoryView {
  isDisabled?: boolean;
  key: CATEGORY_FILTER_ID;
  label: string;
  panels: OntologyPanelCategoryView[];
  selectedViews: SelectCategoryValueView[];
  tooltip?: string;
}

/**
 * View model of single panel view within an ontology-aware filter that contains multiple panels. For example, tissue
 * system, tissue organ or tissue.
 */
export interface OntologyPanelCategoryView {
  label: string;
  isSearchMultiselect: boolean;
  views: SelectCategoryValueView[];
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
 * Model of ontology-aware ancestor/descendant relationships; descandant ontology term IDs keyed by ancestor ontology
 * term ID.
 */
export type OntologyDescendants = { [key: string]: string[] };

/**
 * Prefixes for indicating exact or inferred matches when filtering across category filters that require OR
 * functionality.
 */
export enum OrFilterPrefix {
  EXPLICIT = "E",
  INFERRED = "I",
}

/**
 * Min and max values selected in range category. Empty array if no range is specified (e.g. on clear of range).
 */
export type Range = [number, number] | [];

/**
 * View model of range metadata key.
 */
export interface RangeCategoryView {
  isDisabled?: boolean;
  key: CATEGORY_FILTER_ID;
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
  key: CategoryValueId;
  count: number;
  label: string;
  selected: boolean;
  selectedPartial: boolean;
  value: string; // Value to mark as selected.
}

/**
 * Metadata values grouped by metadata key, for single or multiselect categories.
 */
export interface SelectCategoryView {
  isDisabled?: boolean;
  key: CATEGORY_FILTER_ID;
  label: string;
  pinnedValues: SelectCategoryValueView[];
  tooltip?: string;
  unpinnedValues: SelectCategoryValueView[];
  values: SelectCategoryValueView[]; // both pinned and unpinned values
}

/**
 * Tissue-related filterable values of collections and datasets.
 */
export interface TissueCategories {
  tissueCalculated: string[];
}
