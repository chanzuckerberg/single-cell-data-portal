import { CellValue, FilterValue, Row } from "react-table";
import { EVENTS } from "src/common/analytics/events";
import {
  Collection,
  COLLECTION_STATUS,
  Ontology,
  PublisherMetadata,
} from "src/common/entities";

/**
 * Payload key when tracking select of category values. For example, "organ" in FILTER_SELECT_ORGAN : {organ: "brain"}.
 */
export enum ANALYTICS_PAYLOAD_KEY {
  CELL_CLASS = "cellClass",
  CELL_SUBCLASS = "cellSubclass",
  CELL_TYPE = "cellType",
  ORGAN = "organ",
  SUSPENSION_TYPE = "suspension_type",
  SYSTEM = "system",
  TISSUE = "tissue",
}

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
  isLabelVisible: boolean; // True if ontology label is to be displayed (e.g. "Other Organisms")
  isSearchable: boolean; // True if ontology category search box is displayed
  isZerosVisible: boolean; // True if zero values are to be displayed
  source: OntologyTermSet;
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
  "CELL_TYPE_CALCULATED" = "CELL_TYPE_CALCULATED",
  "CURATOR_NAME" = "CURATOR_NAME",
  "DEVELOPMENT_STAGE" = "DEVELOPMENT_STAGE",
  "DISEASE" = "DISEASE",
  "SELF_REPORTED_ETHNICITY" = "SELF_REPORTED_ETHNICITY",
  "GENE_COUNT" = "GENE_COUNT",
  "ORGANISM" = "ORGANISM",
  "PUBLICATION" = "PUBLICATION",
  "PUBLICATION_DATE_VALUES" = "PUBLICATION_DATE_VALUES",
  "SEX" = "SEX",
  "STATUS" = "STATUS",
  "SUSPENSION_TYPE" = "SUSPENSION_TYPE",
  "TISSUE_CALCULATED" = "TISSUE_CALCULATED",
}

// ** Category filter types ** //

/**
 * Base category filter configuration model containing values shared across the different types of category filters.
 */
export interface BaseCategoryFilterConfig {
  analyticsEvent?: EVENTS;
  analyticsPayloadKey?: ANALYTICS_PAYLOAD_KEY;
  categoryFilterId: CATEGORY_FILTER_ID;
  filterOnKey: FilterKey; // Key in result set row values to filter on.
  label: string;
  multiselect: boolean; // True if category can have multiple values selected.
  pinnedCategoryValues?: CATEGORY_VALUE_KEY[];
  pinnedPosition?: PINNED_POSITION;
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
 */
export type MultiPanelOntologyFilterConfig = BaseCategoryFilterConfig &
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
  | MultiPanelOntologyFilterConfig
  | RangeCategoryFilterConfig
  | SelectCategoryFilterConfig;

// ** Panel search discriminating unions ** //

/**
 * Model of filter that when a search term is selected, the search input is cleared. For example, system or organ.
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
interface CuratedOnlyValueSourceKind {
  source: OntologyTermSet;
  sourceKind: "CURATED";
}

/**
 * Model of a filter panel that displays any value that is explicit only. For example, tissue.
 */
interface ExplicitOnlyValueSourceKind {
  sourceKind: "EXPLICIT_ONLY";
}

// ** Filter value kind discriminating unions ** //

/**
 * Model of a filter panel whose category filter values filter for both inferred and explicit values. For example,
 * tissue system and tissue organ.
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
  analyticsEvent?: EVENTS;
  analyticsPayloadKey?: ANALYTICS_PAYLOAD_KEY;
  id: CATEGORY_FILTER_PANEL_ID;
  label: string;
}

/**
 * Complete set of filter panel IDs.
 */
export enum CATEGORY_FILTER_PANEL_ID {
  "CELL_TYPE_CELL_CLASS" = "CELL_TYPE_CELL_CLASS",
  "CELL_TYPE_CELL_SUBCLASS" = "CELL_TYPE_CELL_SUBCLASS",
  "CELL_TYPE" = "CELL_TYPE",
  "TISSUE_SYSTEM" = "TISSUE_SYSTEM",
  "TISSUE_ORGAN" = "TISSUE_ORGAN",
  "TISSUE" = "TISSUE",
}

/**
 * Filter category panel that uses checkboxes for filtering, is ontology-aware, and whose selected values restrict
 * allowed displayable values in children category filters. For example, tissue system.
 */
type ChildCategoryFilterPanelConfig = BaseCategoryFilterPanelConfig &
  CuratedOnlyValueSourceKind &
  InferredExplicitFilterValueKind &
  SearchSingleSelectKind;

/**
 * Filter category panel that uses checkboxes for filtering, is ontology-aware, and whose selected values restrict
 * allowed displayable values in children category filters. For example, tissue system.
 */
type LeafCategoryFilterPanelConfig = BaseCategoryFilterPanelConfig &
  ExplicitOnlyValueSourceKind &
  ExplicitOnlyFilterValueKind &
  SearchMultiSelect;

/**
 * Filter category panel that uses checkboxes for filtering, is ontology-aware, and whose selected values restrict
 * allowed displayable values in children category filters. For example, tissue system.
 */
type ParentCategoryFilterPanelConfig = BaseCategoryFilterPanelConfig &
  CuratedOnlyValueSourceKind &
  InferredExplicitFilterValueKind &
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
  NO_PUBLICATION = "No publication",
  NORMAL = "normal",
}

/**
 * Filterable values of datasets and collections.
 */
export interface Categories {
  assay: Ontology[];
  cell_count: number;
  cell_type: Ontology[];
  cell_type_ancestors: string[];
  cellTypeCalculated: string[];
  disease: Ontology[];
  development_stage_ancestors: string[];
  self_reported_ethnicity: Ontology[];
  organism: Ontology[];
  sex: Ontology[];
  suspension_type: string[];
  tissue: Ontology[];
  tissue_ancestors: string[];
  tissueCalculated: string[];
}

/**
 * Entry in react-table's filters arrays, models selected category values in a category.
 */
export interface CategoryFilter {
  id: string;
  value: FilterValue;
}

/**
 * Keys of categories that are ontology arrays.
 */
export type CategoriesKeyOfTypeOntologyArray = {
  [K in keyof Categories]: Categories[K] extends Ontology[] ? K : never;
}[keyof Categories];

/*
 * Set of all category values in the full result set, keyed by their corresponding category.
 */
export type CategorySet = { [K in CATEGORY_FILTER_ID]: CategorySetValue };

/**
 * Possible category set values, either a set of category key values (for single or multiselect categories, or ontology
 * categories) or a range.
 */
export type CategorySetValue = Set<CategoryValueId> | Range;

/**
 * Possible category view model types.
 */
export type CategoryView =
  | OntologyCategoryView
  | MultiPanelOntologyCategoryView
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
    Testable {
  curator_name?: string; // Curator view.
  id: string;
  name: string;
  published_at: number;
  publisher_metadata?: PublisherMetadata; // Collection publication metadata
  recency: number; // Used by sort
  revised_at?: number;
  revisedBy?: string; // Curator view.
  status?: COLLECTION_STATUS[]; // Curator view.
  testId: string; // Test ID for e2e tests, facilitates easy look-ups of rows.
}

/**
 * Join of dataset and collection information, optimized for filtering datasets.
 */
export interface DatasetRow
  extends Categories,
    PublisherMetadataCategories,
    Testable {
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
 * Display values of unspecified self-reported ethnicity labels.
 */
export enum SELF_REPORTED_ETHNICITY_UNSPECIFIED_LABEL {
  "unknown" = "Unknown",
}

/**
 * State backing filter functionality and calculations. Converted to view model for display.
 */
export type FilterState = {
  [K in CATEGORY_FILTER_ID]: RangeCategory | KeyedSelectCategoryValue;
};

/**
 * Internal filter model of a single or multiselect category value, or an ontology category value: category value keyed
 * by category value key (for easy look-up when summarizing category).
 */
export type KeyedSelectCategoryValue = Map<
  CategoryValueId,
  SelectCategoryValue
>;

/**
 * Model of filter category configs keyed by category key. Used instead of generic Map to prevent null checking when
 * grabbing keyed value.
 */
export type KeyedCategoryFilterConfigs = {
  [K in CATEGORY_FILTER_ID]: CategoryFilterConfig;
};

/**
 * Model of selected and partially selected states of a multi-panel category filter. Used both internally to record
 * the current selected states, as well as externally when storing multi-panel category filter state to local storage.
 */
export interface MultiPanelCategoryFilterSelectedUIState {
  selected: CategoryValueId[];
  selectedPartial: CategoryValueId[];
}

/**
 * Selected and partially selected states of multi-panel categories filters, keyed by category filter. Returned from
 * hook as an accessor to multi-panel selected state (currently used when saving filter state to local storage).
 */
export type MultiPanelSelectedUIState = {
  [K in CATEGORY_FILTER_ID]?: MultiPanelCategoryFilterSelectedUIState;
};

/**
 * UI model of selected values in a multi-panel category filter. This model contains all selected and partial selected
 * values in a multi-panel category filter and is required as a separate record from react-table's filters which
 * only contains "overridden" selected values. For example, when "digestive system" and "tongue" are both selected in
 * the UI, react-table will only know that "tongue" is selected. It also contains a model of the cross-panel ancestor/
 * descendant relationships.
 */
export interface MultiPanelCategoryFilterUIState
  extends MultiPanelCategoryFilterSelectedUIState {
  uiNodesByCategoryValueId: Map<CategoryValueId, MultiPanelUINode>;
}

/**
 * Model of the cross-panel ancestor/descendant relationships. For example, blood has:
 * - One UI parent (hematopoietic system),
 * - Three UI children (blood non-specific, umbilical cord blood, venous blood).
 * This model facilitates easy traversal and and lookup of ancestor/descendant relationships.
 */
export interface MultiPanelUINode {
  categoryValueId: CategoryValueId;
  panelIndex: number;
  uiChildren: CategoryValueId[];
  uiParents: CategoryValueId[];
}

/**
 * UI model of selected values across all multi-panel category filters.
 */
export type MultiPanelUIState = Map<
  CATEGORY_FILTER_ID,
  MultiPanelCategoryFilterUIState
>;

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
  selectedValue: CategoryValueId | Range,
  selectedLabel: string | Range,
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
  categoryFilterId: CATEGORY_FILTER_ID;
  label: string;
  views: OntologyCategoryTreeView[];
  tooltip?: string;
}

/**
 * View model of an ontology-aware filter that contains multiple panels. For example, tissue.
 */
export interface MultiPanelOntologyCategoryView {
  isDisabled?: boolean;
  categoryFilterId: CATEGORY_FILTER_ID;
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
 * Internal filter model of a range category.
 */
export interface RangeCategory {
  key: CategoryValueId;
  max: number;
  min: number;
  selectedMax?: number;
  selectedMin?: number;
}

/**
 * View model of range metadata key.
 */
export interface RangeCategoryView {
  isDisabled?: boolean;
  categoryFilterId: CATEGORY_FILTER_ID;
  label: string;
  max: number;
  min: number;
  selectedMax?: number;
  selectedMin?: number;
  tooltip?: string;
}

/**
 * Possible locations of pinned values in filter menu.
 */
export enum PINNED_POSITION {
  "BOTTOM" = "BOTTOM",
  "TOP" = "TOP",
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
 * Publication-related filterable values of collections and datasets.
 */
export interface PublisherMetadataCategories {
  // publicationAuthors?: string[];
  publicationDateValues?: number[]; // Set of date bins that publication date falls within
  summaryCitation: string;
}

/**
 * React table count summary values.
 */
export interface TableCountSummary {
  row: number;
  total: number;
}

/**
 * "row" prop passed to react-table's Cell function.
 */
export type RowPropsValue<T extends Categories> = { row: Row<T> };

/**
 * Internal filter model of a single or multiselect category, an ontology category or a multi-panel category.
 */
export interface SelectCategoryValue {
  categoryValueId: CategoryValueId;
  count: number;
  selected: boolean;
  selectedPartial: boolean; // Only applicable to multi-panel categories.
}

/**
 * View model of metadata value, selected state and count for single or multiselect categories.
 */
export interface SelectCategoryValueView {
  categoryValueId: CategoryValueId;
  count: number;
  label: string;
  selected: boolean;
  selectedPartial: boolean;
  visible: boolean;
}

/**
 * Metadata values grouped by metadata key, for single or multiselect categories.
 */
export interface SelectCategoryView {
  isDisabled?: boolean;
  categoryFilterId: CATEGORY_FILTER_ID;
  label: string;
  pinnedPosition?: PINNED_POSITION;
  pinnedValues: SelectCategoryValueView[];
  tooltip?: string;
  unpinnedValues: SelectCategoryValueView[];
  values: SelectCategoryValueView[]; // both pinned and unpinned values
}

/**
 * Interface containing test ID, facilitating easy DOM-lookups of elements in e2e tests.
 */
export interface Testable {
  testId?: string;
}
