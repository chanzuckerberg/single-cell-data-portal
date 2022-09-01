import {
  Dispatch,
  SetStateAction,
  useCallback,
  useEffect,
  useState,
} from "react";
import { Filters, FilterValue, Row } from "react-table";
import { Ontology } from "src/common/entities";
import { TISSUE_DESCENDANTS } from "src/common/queries/tissue-descendants";
import {
  CATEGORY_FILTER_CONFIGS_BY_ID,
  COLLATOR_CASE_INSENSITIVE,
  DEVELOPMENT_STAGE_ONTOLOGY_TERM_SET,
  TISSUE_ORGAN_ONTOLOGY_TERM_SET,
  TISSUE_SYSTEM_ONTOLOGY_TERM_SET,
} from "src/components/common/Filter/common/constants";
import {
  Categories,
  CategoryFilterConfig,
  CategoryFilterPanelConfig,
  CategoryValueId,
  CategoryView,
  CATEGORY_FILTER_ID,
  CuratedOntologyCategoryFilterConfig,
  ETHNICITY_UNSPECIFIED_LABEL,
  FilterKey,
  OnFilterFn,
  OntologyCategoryTreeNodeView,
  OntologyCategoryTreeView,
  OntologyCategoryView,
  OntologyMultiPanelCategoryView,
  OntologyMultiPanelFilterConfig,
  OntologyNode,
  OntologyPanelCategoryView,
  OntologyTermSet,
  ONTOLOGY_VIEW_KEY,
  ONTOLOGY_VIEW_LABEL,
  ON_FILTER_SOURCE,
  OrFilterPrefix,
  ORGANISM,
  PUBLICATION_DATE_LABELS,
  Range,
  RangeCategoryView,
  SelectCategoryValueView,
  SelectCategoryView,
} from "src/components/common/Filter/common/entities";
import {
  buildExplicitOntologyTermId,
  buildInferredOntologyTermId,
  findOntologyNodeById,
  findOntologyParentNode,
  getOntologySpeciesKey,
  listOntologyTreeIds,
  removeOntologyTermId,
  removeOntologyTermIdPrefix,
  splitOntologyTermIdAndPrefix,
} from "src/components/common/Filter/common/utils";
import { track } from "../analytics";

// TODO(cc) separate to utils files.
// TODO(cc) publication date sort order

/**
 * String value to append to labels in multi-panel categories if the value appears in more than one panel.
 */
const LABEL_SUFFIX_NON_SPECIFIC = ", non-specific";

/**
 * Entry in react-table's filters arrays, models selected category values in a category.
 */
interface CategoryFilter {
  id: string;
  value: FilterValue;
}

/*
 * Set of all category values in the full result set, keyed by their corresponding category.
 */
type CategorySet = { [K in CATEGORY_FILTER_ID]: CategorySetValue };

/**
 * Possible category set values, either a set of category key values (for single or multiselect categories, or ontology
 * categories) or a range.
 */
type CategorySetValue = Set<CategoryValueId> | Range;

/**
 * Internal filter model of a single or multiselect category value, or an ontology category value: category value keyed
 * by category value key (for easy look-up when summarizing category).
 */
type KeyedSelectCategoryValue = Map<CategoryValueId, SelectCategoryValue>;

/**
 * Internal filter model of a single or multiselect category, an ontology category or a multi-panel category.
 */
interface SelectCategoryValue {
  key: CategoryValueId;
  count: number;
  selected: boolean;
  selectedPartial: boolean; // Only applicable to multi-panel categories.
}

/**
 * Shape of return value from this useFilter hook.
 */
export interface FilterInstance {
  categoryViews: CategoryView[];
  onFilter: OnFilterFn;
}

/**
 * State backing filter functionality and calculations. Converted to view model for display.
 */
type FilterState = {
  [K in CATEGORY_FILTER_ID]: RangeCategory | KeyedSelectCategoryValue;
};

/**
 * UI model of selected values in a multi-panel category filter. This model contains all selected and partial selected
 * values in a multi-panel category filter and is required as a separate record from react-table's filters which
 * only contains "overridden" selected values. For example, when "digestive system" and "tongue" are both selected in
 * the UI, react-table will only know that "tongue" is selected.
 * TODO(cc) rename, docs
 */
export interface MultiPanelCategoryFilterUIState {
  selected: CategoryValueId[];
  selectedPartial: CategoryValueId[];
  uiNodesByCategoryValueId: Map<CategoryValueId, MultiPanelUINode>;
}

// TODO(cc) - move to, export here and throughout
export interface MultiPanelUINode {
  categoryValueId: CategoryValueId;
  uiChildren: CategoryValueId[];
  uiParents: CategoryValueId[];
}

/**
 * UI model of selected values across all multi-panel category filters.
 */
type MultiPanelUIState = Map<
  CATEGORY_FILTER_ID,
  MultiPanelCategoryFilterUIState
>;

/**
 * Selected filters applicable to a category; used when deriving category value counts from current set of filters.
 * Identical queries can be shared by categories to reduce the number of result set filtering.
 */
interface Query<T extends Categories> {
  categoryFilterIds: CATEGORY_FILTER_ID[];
  filters: Filters<T>;
}

/**
 * Internal filter model of a range category.
 */
interface RangeCategory {
  key: CategoryValueId;
  max: number;
  min: number;
  selectedMax?: number;
  selectedMin?: number;
}

/**
 * react-table function to call when updating set of selected filters.
 */
/* eslint-disable-next-line @typescript-eslint/no-explicit-any -- function type as per react-table's setFilter. */
type SetFilterFn = (columnId: string, updater: any) => void;

/**
 * Generic tooltip displayed when select, ontology or range category is disabled due to no values matching current
 * filter.
 */
const TOOLTIP_CATEGORY_DISABLED =
  "There are no values that meet the current filter criteria";

/**
 * Faceted filter functionality over dataset metadata. "or" between values, "and" across categories.
 * @param originalRows - Original result set before filtering.
 * @param categoryFilterIds - Set of category IDs to include for this filter instance.
 * @param filters - Current set of selected category values (values) keyed by category (id).
 * @param setFilter - Function to update set of selected values for a category.
 * @returns Object containing filter accessor (view model of filter state) and filter mutator (function to modify react-
 * table's internal filter state).
 */
export function useCategoryFilter<T extends Categories>(
  originalRows: Row<T>[],
  categoryFilterIds: Set<CATEGORY_FILTER_ID>,
  filters: Filters<T>,
  setFilter: SetFilterFn
): FilterInstance {
  // Complete set of categories and category values for the result set.
  const [categorySet, setCategorySet] = useState<CategorySet>();

  // Core filter state facilitating build of complete set of categories, category values and counts for a filtered
  // result set.
  const [filterState, setFilterState] = useState<FilterState>();

  // Set of ontology term labels keyed by ontology term ID required for label lookup when building view models for
  // ontology term-backed fields.
  const [ontologyTermLabelsById, setOntologyTermLabelsById] =
    useState<Map<string, string>>();

  // Internally saved selected values for each multi-panel category filter; used to set selected state of values in
  // ontology-aware category filters. This is required for category filters where cross-panel restrictions are applied
  // (e.g. tissue system restricts tissue organ and tissue). We can not use react-table's filters as it only contains
  // the most restrictive value (e.g. if "renal system" and "kidney" are both selected, only "kidney" is set as a
  // selected value in react-table). We need a variable to save *all* selected values so this can be reflected in the
  // view models. Also saved is the delta between the multi-panel selected values and the selected values passed to
  // react-table, facilitating easy calculation of partially selected values and cross-category view restrictions.
  const [multiPanelUIState, setMultiPanelUIState] = useState<
    Map<CATEGORY_FILTER_ID, MultiPanelCategoryFilterUIState>
  >(new Map());

  // Set up original, full set of categories and their values.
  useEffect(() => {
    // Only build category set if there are rows to parse category values from. Only build category set once on load.
    if (!originalRows.length || categorySet) {
      return;
    }

    setCategorySet(buildCategorySet(originalRows, categoryFilterIds));
  }, [originalRows, categoryFilterIds, categorySet]);

  // Build up map of ontology term labels keyed by ID.
  useEffect(() => {
    // Only build category set if there are rows to parse category values from. Only build category set once on load.
    if (!originalRows.length || ontologyTermLabelsById) {
      return;
    }

    setOntologyTermLabelsById(keyOntologyTermLabelsById(originalRows));
  }, [originalRows, categoryFilterIds, ontologyTermLabelsById]);

  // Build up UI hierarchies for each multi-panel category filter, used to facilitate easy calculation of selected
  // and partially selected states.
  useEffect(() => {
    // Only set multi-panel state if there are rows to parse category values from. TODO(cc) check only called once
    if (!originalRows.length) {
      return;
    }

    setMultiPanelUIState(
      buildMultiPanelUIState(originalRows, categoryFilterIds)
    );
  }, [categoryFilterIds, originalRows]);

  // Build next filter state on change of filter.
  useEffect(() => {
    // Must have category set before next filter state can be calculated.
    if (!categorySet) {
      return;
    }

    const nextFilterState = buildNextFilterState(
      originalRows,
      categoryFilterIds,
      filters,
      multiPanelUIState,
      categorySet
    );
    setFilterState(nextFilterState);
  }, [
    categoryFilterIds,
    categorySet,
    filters,
    originalRows,
    multiPanelUIState,
  ]);

  // Update set of filters on select of category value. Track selected category value.
  const onFilter = useCallback<OnFilterFn>(
    (
      categoryFilterId: CATEGORY_FILTER_ID,
      categoryValueKey: CategoryValueId | null, // TODO(cc) is this still necessary
      selectedValue: CategoryValueId | Range,
      source: ON_FILTER_SOURCE = ON_FILTER_SOURCE.FILTER
    ) => {
      if (!categorySet) {
        return; // Error state - category set should be set at this point.
      }

      // Grab the configuration model for the selected category.
      const config = CATEGORY_FILTER_CONFIGS_BY_ID[categoryFilterId];

      // Handle range categories.
      if (isSelectedValueRange(selectedValue)) {
        onFilterRangeCategory(config, selectedValue, setFilter);
        return;
      }

      if (!categoryValueKey) {
        return; // Error state - category value key must be defined for select and ontology categories.
      }

      // Handle ontology categories.
      if (isCuratedOntologyCategoryFilterConfig(config)) {
        onFilterCuratedOntologyCategory(
          config,
          categoryValueKey,
          selectedValue,
          filters,
          setFilter,
          categorySet
        );
        return;
      }

      // Handle multi-panel categories. TODO(cc) don't pass in setFilter, here and for others.
      if (isMultiPanelCategoryFilterConfig(config)) {
        onFilterMultiPanelCategory(
          config,
          categoryValueKey,
          selectedValue,
          source,
          setFilter,
          multiPanelUIState,
          setMultiPanelUIState
        );
        return;
      }

      // Handle single or multiselect categories.
      onFilterSelectCategory(
        config,
        categoryValueKey,
        selectedValue,
        filters,
        setFilter
      );
    },
    [categorySet, filters, setFilter, multiPanelUIState]
  );

  return {
    categoryViews: buildCategoryViews(
      filterState,
      multiPanelUIState,
      ontologyTermLabelsById
    ),
    onFilter,
  };
}

/**
 * Add back any category values that have been filtered out, set their values to 0.
 * @param nextFilterState - Filter state currently being built due to change in filter.
 * @param categorySet - Original, unfiltered sets of category values keyed by their category.
 */
function addEmptyCategoryValues(
  nextFilterState: FilterState,
  categorySet: CategorySet
) {
  // Check filter state for each category for missing (empty) category values.
  for (const [categoryFilterId, categoryValuesByKey] of Object.entries(
    nextFilterState
  )) {
    // Adding back empty category values is only applicable to select or ontology category values.
    if (!isSelectCategoryValue(categoryValuesByKey)) {
      continue;
    }

    // Grab the expected set of category values.
    const allCategoryValueKeys =
      categorySet[categoryFilterId as CATEGORY_FILTER_ID];
    if (
      !allCategoryValueKeys || // Error state - all category values for this category can't be found.
      !isCategoryValueIdSet(allCategoryValueKeys) // Error state - should be category key value.
    ) {
      return;
    }

    // If expected category value is missing from this category's category values, add it back in with a count of 0.
    [...allCategoryValueKeys.values()].forEach(
      (categoryValueKey: CategoryValueId) => {
        if (!categoryValuesByKey.has(categoryValueKey)) {
          categoryValuesByKey.set(categoryValueKey, {
            count: 0,
            key: categoryValueKey,
            selected: false,
            selectedPartial: false,
          });
        }
      }
    );
  }
}

/**
 * Add all descendents to the set of selected values. Do not add IDs of descendents that are specified in the
 * ontology view but not present in the full set of data as this will filter the data by values that do not exist
 * (resulting in an empty set).
 * @param ontologyNodes - Nodes to add to selected values.
 * @param selectedCategoryValues - The current set of selected values.
 * @param categoryKeyValues - Original, full set of values for this category.
 */
function addOntologyDescendents(
  ontologyNodes: OntologyNode[],
  selectedCategoryValues: Set<CategoryValueId>,
  categoryKeyValues: Set<CategoryValueId>
) {
  // For each node, add self and possibly children.
  ontologyNodes.forEach((ontologyNode) => {
    const { ontology_term_id: ontologyId } = ontologyNode;
    // Only add self if memeber of original set of values for this category.
    if (categoryKeyValues.has(ontologyId)) {
      selectedCategoryValues.add(ontologyId);
    }
    // Add children, if any.
    if (ontologyNode.children) {
      addOntologyDescendents(
        ontologyNode.children,
        selectedCategoryValues,
        categoryKeyValues
      );
    }
  });
}

/**
 * Add range categories to the next filter state. Range categories are not summarized and must be explicitly added from
 * the original category set state.
 * @param categoryFilterIds - Set of category IDs to include for this filter instance.
 * @param nextFilterState - Filter state currently being built due to change in filter.
 * @param categorySet - Original, unfiltered sets of category values keyed by their category.
 */
function addRangeCategories(
  categoryFilterIds: Set<CATEGORY_FILTER_ID>,
  nextFilterState: FilterState,
  categorySet: CategorySet
) {
  Array.from(categoryFilterIds.values()).forEach((categoryKey) => {
    // Grab the expected range for this category.
    const categorySetRange = categorySet[categoryKey];
    if (
      !categorySetRange || // Error state - original range for this category can't be found.
      isCategoryValueIdSet(categorySetRange) // Error state - should be a range.
    ) {
      return;
    }

    // Add range to next filter state.
    const [min, max] = categorySetRange;
    nextFilterState[categoryKey] = {
      key: categoryKey,
      max: max ?? 0,
      min: min ?? 0,
    };
  });
}

/**
 * Determine the filter type for each selected filter value and apply to rows.
 * @param originalRows - Original result set before filtering.
 * @param filters - Selected filters to apply to rows.
 * @returns Filtered array of rows.
 */
function applyFilters<T extends Categories>(
  originalRows: Row<T>[],
  filters: Filters<T>
): Row<T>[] {
  // Return all rows if there are no filters.
  if (filters.length === 0) {
    return originalRows;
  }
  return originalRows.filter((row: Row<T>) => {
    // "and" across categories.
    return filters.every((filter: CategoryFilter) => {
      const rowValue = row.values[filter.id];
      if (isMatchKindBetween(filter.id as CATEGORY_FILTER_ID)) {
        return between(rowValue, filter);
      }
      return includesSome(rowValue, filter);
    });
  });
}

/**
 * Determine the rows that have values between the given filter value. Mimics react-query's between functionality.
 * @param rowValue - Value to filter row by.
 * @param filter - Selected filter to apply to row.
 * @returns Filtered array of rows.
 */
function between(rowValue: string | string[], filter: CategoryFilter): boolean {
  const [min, max] = filter.value;
  return Boolean(rowValue && rowValue >= min && rowValue <= max);
}

/**
 * Set up model of original, complete set of categories and their values.
 * @param originalRows - Original result set before filtering.
 * @param categoryFilterIds - Set of category filter IDs to include for this filter instance.
 * @returns Sets of category values keyed by their category.
 */
function buildCategorySet<T extends Categories>(
  originalRows: Row<T>[],
  categoryFilterIds: Set<CATEGORY_FILTER_ID>
): CategorySet {
  // Build up category values for each category
  return Array.from(categoryFilterIds.values()).reduce(
    (accum: CategorySet, categoryFilterId: CATEGORY_FILTER_ID) => {
      // Calculate the initial state of range categories.
      if (isMatchKindBetween(categoryFilterId)) {
        const counts = originalRows
          // Use filterKey to pull value from original row
          .map((originalRow) => originalRow.values[categoryFilterId])
          .filter((count) => !!count || count === 0); // Remove bad data, just in case!

        accum[categoryFilterId] = [
          counts.length ? Math.min(...counts) : 0,
          counts.length ? Math.max(...counts) : 0,
        ];
        return accum;
      }

      // Determine the set of ontology IDs in the ontology tree for this category, if applicable. There are possibly
      // more ontology IDs listed for each row than we want to display (that is, a row can possibly have a higher
      // granularity of ontology IDs than the UI is to display). We'll use the ontology IDs of the ontology tree
      // to determine which row values are to be included in the category set.
      const config = CATEGORY_FILTER_CONFIGS_BY_ID[categoryFilterId];
      const isCategoryOntology = isCuratedOntologyCategoryFilterConfig(config);
      let categoryOntologyIds: Set<string>;
      if (isCategoryOntology) {
        categoryOntologyIds = listOntologyTreeIds(config.mask);
      }

      // Handle single or multi select categories. Check category value for this category, in every row.
      originalRows.forEach((originalRow: Row<T>) => {
        // Grab the category values already added for this category, create new set if it hasn't already been created.
        let categoryValueSet = accum[categoryFilterId] as Set<CategoryValueId>;
        if (!categoryValueSet) {
          categoryValueSet = new Set<CategoryValueId>();
          accum[categoryFilterId] = categoryValueSet;
        }
        // Add the category values for this row to the set.
        let values: CategoryValueId | CategoryValueId[] =
          originalRow.values[categoryFilterId];
        if (typeof values === "undefined") {
          console.log(`No values found for category "${categoryFilterId}".`);
          return accum;
        }
        if (!Array.isArray(values)) {
          values = [values];
        }

        // Add each value to category set.
        values
          .filter((value: CategoryValueId) => {
            // If category is an ontology, confirm value is included in ontology tree for display.
            return !isCategoryOntology || categoryOntologyIds.has(value);
          })
          .forEach((value: CategoryValueId) => categoryValueSet.add(value));
      });
      return accum;
    },
    {} as CategorySet
  );
}

/**
 * Build the display value for the given category and category value. For ontology terms, look up corresponding labels.
 * @param config - Config model of category to build category value views for.
 * @param categoryValueKey - Category value to display (e.g. "normal").
 * @param ontologyTermLabelsById - Set of ontology term labels keyed by term ID, used to determine labels for ontology
 * @returns String to display as a label for the given category and category value.
 */
function buildCategoryValueLabel(
  config: CategoryFilterConfig,
  categoryValueKey: CategoryValueId,
  ontologyTermLabelsById: Map<string, string>
): string {
  const { categoryFilterId } = config;
  if (categoryFilterId === CATEGORY_FILTER_ID.PUBLICATION_DATE_VALUES) {
    return PUBLICATION_DATE_LABELS[
      `LABEL_${categoryValueKey}` as keyof typeof PUBLICATION_DATE_LABELS
    ];
  }

  if (
    categoryFilterId === CATEGORY_FILTER_ID.ETHNICITY &&
    !isEthnicitySpecified(categoryValueKey)
  ) {
    return ETHNICITY_UNSPECIFIED_LABEL[
      categoryValueKey as keyof typeof ETHNICITY_UNSPECIFIED_LABEL
    ];
  }

  // Look up labels for ontology term IDs.
  if (config.labelKind === "LOOKUP_LABEL_BY_TERM_ID") {
    let processedCategoryValueKey = categoryValueKey;
    if (config.categoryFilterId === "TISSUE_CALCULATED") {
      processedCategoryValueKey = removeOntologyTermIdPrefix(categoryValueKey);
    }

    return (
      ontologyTermLabelsById.get(processedCategoryValueKey) ?? categoryValueKey
    );
  }

  // Return all other category values as is.
  return categoryValueKey;
}

/**
 * Build view-specific models from filter state, to facilitate easy rendering.
 * @param filterState - Categories, category value and their counts with the current filter applied.
 * @param multiPanelUIState - Current set of category values that the user has selected in multi-panel cateogry filters.
 * @param ontologyTermLabelsById - Set of ontology term labels keyed by term ID, used to determine labels for ontology
 * terms.
 * @returns Array of category views objects.
 */
function buildCategoryViews(
  filterState?: FilterState,
  multiPanelUIState?: MultiPanelUIState,
  ontologyTermLabelsById?: Map<string, string>
): CategoryView[] {
  if (!filterState || !multiPanelUIState || !ontologyTermLabelsById) {
    return [];
  }
  return Object.keys(filterState)
    .map((strCategoryFilterId: string) => {
      const categoryFilterId = strCategoryFilterId as CATEGORY_FILTER_ID;

      // Return if there's no filter state for the UI config. This will be true for categories such as cell count and
      // gene count when viewing collections.
      const rangeOrSelectValue = filterState[categoryFilterId];

      // Grab the config for this category filter.
      const config = CATEGORY_FILTER_CONFIGS_BY_ID[categoryFilterId];

      // Build view model for single or multiselect categories
      // TODO(cc) remove need for multiple checks, can we get just kind to assert the category value type?
      if (
        config.viewKind === "SELECT" &&
        isSelectCategoryValue(rangeOrSelectValue)
      ) {
        return buildSelectCategoryView(
          config,
          rangeOrSelectValue,
          filterState,
          ontologyTermLabelsById
        );
      }

      // Build view model for curated ontology categories.
      if (
        config.viewKind === "CURATED_ONTOLOGY" &&
        isSelectCategoryValue(rangeOrSelectValue) &&
        isCuratedOntologyCategoryFilterConfig(config)
      ) {
        return buildCuratedOntologyCategoryView(
          config,
          rangeOrSelectValue,
          filterState,
          ontologyTermLabelsById
        );
      }

      // Build view model for range categories.
      if (
        config.viewKind === "RANGE" &&
        !isSelectCategoryValue(rangeOrSelectValue)
      ) {
        return buildRangeCategoryView(config, rangeOrSelectValue);
      }

      // Build view model for multi-panel categories.
      if (
        config.viewKind === "MULTI_PANEL" &&
        isSelectCategoryValue(rangeOrSelectValue)
      ) {
        return buildMultiPanelCategoryView(
          config,
          rangeOrSelectValue,
          multiPanelUIState,
          ontologyTermLabelsById
        );
      }

      // Error - unknown view kind/category value.
      console.log(
        `Error attempting to build view model for category "${config.categoryFilterId}".`
      );
    })
    .filter((categoryView): categoryView is CategoryView => !!categoryView)
    .sort(sortCategoryViews);
}

/**
 * TODO(cc) docs, location, separate functions
 */
function buildMultiPanelUIState<T extends Categories>(
  originalRows: Row<T>[],
  categoryFilterIds: Set<CATEGORY_FILTER_ID>
): MultiPanelUIState {
  return Array.from(categoryFilterIds.values()).reduce(
    (
      accum: Map<CATEGORY_FILTER_ID, MultiPanelCategoryFilterUIState>,
      categoryFilterId: CATEGORY_FILTER_ID
    ) => {
      // Ignore categories that are not multi-panels.
      const config = CATEGORY_FILTER_CONFIGS_BY_ID[categoryFilterId];
      if (!isMultiPanelCategoryFilterConfig(config)) {
        return accum;
      }

      // Determine the set of values for each panel.
      const categoryValueIdsByPanel = keyCategoryValueIdsByPanel(
        config,
        originalRows
      );

      // Iterate over each set of panel values and build up parents and children for each value.
      const uiHierarchyByCategoryValue = buildUINodesByCategoryValueId(
        categoryValueIdsByPanel
      );

      accum.set(categoryFilterId, {
        selected: [],
        selectedPartial: [],
        uiNodesByCategoryValueId: uiHierarchyByCategoryValue,
      });

      return accum;
    },
    new Map<CATEGORY_FILTER_ID, MultiPanelCategoryFilterUIState>()
  );
}

/**
 * TODO(cc)
 */
export function buildUINodesByCategoryValueId(
  categoryValueIdsByPanel: CategoryValueId[][]
): Map<CategoryValueId, MultiPanelUINode> {
  return categoryValueIdsByPanel.reduce(
    (
      uiAccum: Map<CategoryValueId, MultiPanelUINode>,
      categoryValueIds: CategoryValueId[],
      index: number
    ) => {
      // Add all values in this panel to the map of values.
      categoryValueIds.forEach((categoryValueId) =>
        uiAccum.set(categoryValueId, {
          categoryValueId: categoryValueId,
          uiChildren: [],
          uiParents: [],
        })
      );

      // Build parent child relationships for each value in this panel.
      const parentCategoryValueIdsByPanel = categoryValueIdsByPanel.slice(
        0,
        index
      );
      categoryValueIds.forEach((categoryValueId) => {
        linkParentsAndChildren(
          categoryValueId,
          parentCategoryValueIdsByPanel,
          uiAccum
        );
      });

      return uiAccum;
    },
    new Map<CategoryValueId, MultiPanelUINode>()
  );
}

/**
 * TODO(cc)
 */
export function keyCategoryValueIdsByPanel<T extends Categories>(
  config: OntologyMultiPanelFilterConfig,
  originalRows: Row<T>[]
): CategoryValueId[][] {
  const { panels: panelConfigs } = config;
  return panelConfigs.reduce(
    (uiAccum: CategoryValueId[][], panelConfig: CategoryFilterPanelConfig) => {
      let prefixedOntologyTermIds;

      // Determine the set of values for curated ontology panels.
      if (panelConfig.sourceKind === "CURATED_CATEGORIES") {
        prefixedOntologyTermIds = [
          ...listOntologyTreeIds(panelConfig.mask),
        ].map((ontologyTermId) => buildInferredOntologyTermId(ontologyTermId));
      }
      // Otherwise, build up the set of values for this panel from the original rows.
      else {
        // TODO(cc) reduce rather than create and set/array
        prefixedOntologyTermIds = [
          ...new Set(
            originalRows
              // eslint-disable-next-line @typescript-eslint/ban-ts-comment --- TODO(cc) revisit - different between FilterKey (any value from dataset or collection) vs T extends Categories, also see if type assertion can be resolved
              // @ts-ignore
              .map((originalRow) => originalRow.original[config.filterOnKey]) // TODO(cc)
              .flat()
              .filter(
                (value) =>
                  removeOntologyTermId(value) === OrFilterPrefix.EXPLICIT
              )
          ),
        ];
      }

      uiAccum.push(prefixedOntologyTermIds);

      return uiAccum;
    },
    [] as CategoryValueId[][]
  );
}

/**
 * TODO(cc)
 */
function linkParentsAndChildren(
  categoryValueId: CategoryValueId,
  parentCategoryValueIdsByPanel: CategoryValueId[][],
  uiNodesByCategoryValueId: Map<CategoryValueId, MultiPanelUINode>
) {
  // Link parent and children between this panel and the parent panel if:
  // 1. Child value is a descendant of the parent value and,
  // 2. Child value is not a descendant of any children already added to parent.
  parentCategoryValueIdsByPanel.forEach((parentCategoryValueIds) => {
    parentCategoryValueIds.forEach((panelCategoryValueId) => {
      // Check if value is a descendant of panel value.
      const isDescendantOfPanelValue = isDescendant(
        categoryValueId,
        panelCategoryValueId,
        TISSUE_DESCENDANTS
      );

      // If value isn't a descendant, there's no parent child relationship to link here.
      if (!isDescendantOfPanelValue) {
        return;
      }

      // Value is a descendant, check it is not "blocked". That is, the value is not a descendant of any children of
      // the panel value.
      // Note: the logic here possibly might need revisiting if more than three panels are added. For example, is it
      // possible that the children list is ever incomplete at this point?
      // TODO(cc) revisit chaining here and below
      const panelCategoryValueUIChildren =
        uiNodesByCategoryValueId.get(panelCategoryValueId)?.uiChildren ?? [];
      const isDescendantOfUIChild = panelCategoryValueUIChildren.some(
        (panelCategoryValueUIChild) => {
          // There are no descendants of explicit values; value can't be a descendant if child is explicit.
          if (isExplicitOntologyTermId(panelCategoryValueUIChild)) {
            return false;
          }

          // Check if value is a descendant of panel children.
          return isDescendant(
            categoryValueId,
            panelCategoryValueUIChild,
            TISSUE_DESCENDANTS
          );
        }
      );

      // Only add value if it's not "blocked".
      if (!isDescendantOfUIChild) {
        linkParentAndChild(
          categoryValueId,
          panelCategoryValueId,
          uiNodesByCategoryValueId
        );
      }
    });
  });
}

/**
 * TODO(cc) optional chaining
 */
function linkParentAndChild(
  categoryValueId: CategoryValueId,
  parentCategoryValueId: CategoryValueId,
  uiNodesByCategoryValueId: Map<CategoryValueId, MultiPanelUINode>
) {
  uiNodesByCategoryValueId
    .get(parentCategoryValueId)
    ?.uiChildren?.push(categoryValueId);
  uiNodesByCategoryValueId
    ?.get(categoryValueId)
    ?.uiParents.push(parentCategoryValueId);
}

/**
 * TODO(cc) check strings here, check usages of TISSUE_DESCENDANTS (can we put the stripping of the prefix in here?), type for descendants here and below
 */
function isDescendant(
  categoryValueId: CategoryValueId,
  panelCategoryValueId: CategoryValueId,
  descendants: { [p: string]: string[] }
): boolean {
  const ontologyTermId = removeOntologyTermIdPrefix(categoryValueId);
  const panelOntologyTermId = removeOntologyTermIdPrefix(panelCategoryValueId);
  const isDescendantTerm = (descendants[panelOntologyTermId] ?? []).includes(
    ontologyTermId
  );

  const isExplicitTerm = isExplicitTermOfInferredTerm(
    categoryValueId,
    panelCategoryValueId
  );

  return isDescendantTerm || isExplicitTerm;
}

/**
 * TODO(cc) move to utils?
 */
function isExplicitTermOfInferredTerm(
  categoryValueId: CategoryValueId,
  panelCategoryValueId: CategoryValueId
): boolean {
  const ontologyTermId = removeOntologyTermIdPrefix(categoryValueId);
  const panelOntologyTermId = removeOntologyTermIdPrefix(panelCategoryValueId);
  return (
    isExplicitOntologyTermId(categoryValueId) &&
    !isExplicitOntologyTermId(panelCategoryValueId) &&
    ontologyTermId === panelOntologyTermId
  );
}

/**
 * TODO(cc) move to utils?
 */
function isExplicitOntologyTermId(categoryValueId: CategoryValueId): boolean {
  const prefix = removeOntologyTermId(categoryValueId);
  return prefix === OrFilterPrefix.EXPLICIT;
}

/**
 * Build categories, category values and counts for the updated filter. For each category, build up category values
 * counts by counting occurrences of category values across rows. Maintain selected category values state from filters.
 * Retain category values with 0 counts from given category set.
 * @param originalRows - Original result set before filtering.
 * @param categoryFilterIds - Set of category IDs to include for this filter instance.
 * @param filters - Current set of selected category values (values) keyed by category (id).
 * @param multiPanelUIState - Current set of category values that the user has selected on the UI for all multi-
 * panel category filters.
 * @param categorySet - Original, unfiltered sets of category values keyed by their category.
 * @returns New filter state generated from the current set of selected category values.
 */
function buildNextFilterState<T extends Categories>(
  originalRows: Row<T>[],
  categoryFilterIds: Set<CATEGORY_FILTER_ID>,
  filters: Filters<T>,
  multiPanelUIState: MultiPanelUIState,
  categorySet: CategorySet
): FilterState {
  // Build set of filters that are applicable to each category.
  const queries = buildQueries(categoryFilterIds, filters);

  // Build up base filter state of categories, category values and counts.
  const nextFilterState = summarizeCategories(originalRows, queries);

  // Always display category values even if their count is 0; add back any category values that have been filtered out.
  // This allows us to maintain the selected state of values that are selected but possibly filtered out.
  addEmptyCategoryValues(nextFilterState, categorySet);

  // Add range categories to next filter state.
  addRangeCategories(categoryFilterIds, nextFilterState, categorySet);

  // Update selected flag for the selected category values, or selected ranged for range categories.
  setSelectedStates(nextFilterState, filters, multiPanelUIState);

  return nextFilterState;
}

/**
 * Build updated set of selected filters for the given ontology tree category and the selected category value.
 * @param categoryFilterId - Category ID (i.e. "development stage") of selected category value.
 * @param selectedValue - Selected category value ID (e.g. "HsapDv:0000003") to update selected state of.
 * @param filters - Current set of selected category values.
 * @param categoryKeyValues - Original, full set of values for this category.
 * @param mask - View model of ontology for this category.
 * @returns Array of selected category values for the given category.
 */
export function buildNextOntologyCategoryFilters<T extends Categories>(
  categoryFilterId: CATEGORY_FILTER_ID,
  selectedValue: CategoryValueId,
  filters: Filters<T>,
  categoryKeyValues: Set<CategoryValueId>,
  mask: OntologyTermSet
): CategoryValueId[] {
  // Grab the current selected values for the category.
  const categoryFilters = new Set(
    getCategoryFilter(categoryFilterId, filters)?.value as CategoryValueId[]
  );

  // Find the selected and parent node, if any, for the selected value.
  const ontologySpeciesKey = getOntologySpeciesKey(selectedValue);
  const ontologyRootNodes = mask[ontologySpeciesKey];
  if (!ontologyRootNodes) {
    return [...categoryFilters.values()]; // Error state - ontology does not exist.
  }
  const selectedOntologyNode = findOntologyNodeById(
    ontologyRootNodes,
    selectedValue
  );
  if (!selectedOntologyNode) {
    return [...categoryFilters.values()]; // Error state - ontology node with given ID does not exist.
  }
  const parentNode = findOntologyParentNode(
    ontologyRootNodes,
    selectedOntologyNode
  );

  // Toggle selected state of selected category value.
  if (categoryFilters.has(selectedValue)) {
    // Selected value is already in the set of selected values, remove it.
    categoryFilters.delete(selectedValue);

    // Also remove any descendents from the selected set.
    if (selectedOntologyNode.children) {
      removeOntologyDescendents(selectedOntologyNode.children, categoryFilters);
    }

    // Reevaluate selected state of parent. If all children were previously selected and now only some are selected,
    // then parent should no longer be selected.
    if (parentNode) {
      handleOntologyChildRemoved(
        ontologyRootNodes,
        parentNode,
        categoryFilters
      );
    }
  } else {
    // Add selected value to selected set.
    categoryFilters.add(selectedValue);

    // Add all descendents of selected value, if any.
    if (selectedOntologyNode.children) {
      addOntologyDescendents(
        selectedOntologyNode.children,
        categoryFilters,
        categoryKeyValues
      );
    }

    // Reevaluate selected state of parent. If all children are now selected, then parent should also be selected.
    if (parentNode) {
      handleOntologyChildAdded(
        ontologyRootNodes,
        parentNode,
        categoryFilters,
        categoryKeyValues
      );
    }
  }

  return [...categoryFilters.values()];
}

/**
 * Build updated set of selected filters for the given single or multiselect category and the selected category values.
 * @param config - Configuration model of selected category.
 * @param selectedValue - Category value ID to toggle the selected state of.
 * @param filters - Current set of selected category values.
 * @returns Array of selected category values for the given category.
 */
function buildNextSelectCategoryFilters<T extends Categories>(
  config: CategoryFilterConfig,
  selectedValue: CategoryValueId,
  filters: Filters<T>
): CategoryValueId[] {
  const { categoryFilterId, multiselect } = config;

  // Grab the current selected values for the category.
  const categoryFilters = getCategoryFilter(categoryFilterId, filters);

  // Currently, no filters already selected for this category; add category value as first.
  if (!categoryFilters) {
    return [selectedValue];
  }

  // Create new array of selected category value keys, with the selected state of the given category value toggled.
  return toggleCategoryValueSelected(
    selectedValue,
    categoryFilters.value,
    multiselect
  );
}

/**
 * Determine the set of filters that are applicable to each category. That is, for a category, all selected filters
 * other than the selected filters for that category can be applied to the result set to determine the counts for.
 * @param categoryFilterIds - Set of category IDs to include for this filter instance.
 * @param filters - Current set of selected category values (values) keyed by category (id).
 * @returns Array of query models representing of the selected filters applicable for each category.
 */
function buildQueries<T extends Categories>(
  categoryFilterIds: Set<CATEGORY_FILTER_ID>,
  filters: Filters<T>
): Query<T>[] {
  return Array.from(categoryFilterIds.values()).reduce(
    (accum: Query<T>[], categoryFilterId: CATEGORY_FILTER_ID) => {
      // Determine the filters that are applicable to this category.
      const applicableFilters = filters.filter(
        (filter: CategoryFilter) => filter.id !== categoryFilterId
      );

      // Check if we have an existing query with an identical filter. If so, add category to that query. Otherwise,
      // create new query for this filter.
      const matchingQuery = accum.find((query: Query<T>) =>
        isFilterEqual(query.filters, applicableFilters)
      );
      if (matchingQuery) {
        (matchingQuery as Query<T>).categoryFilterIds.push(categoryFilterId);
      } else {
        accum.push({
          categoryFilterIds: [categoryFilterId],
          filters: applicableFilters,
        });
      }
      return accum;
    },
    []
  );
}

/**
 * Build view model of a curated ontology category such as development stage.
 * @param config - Config model of a curated ontology category.
 * @param categoryValueByValue - Internal filter model of ontology category.
 * @param filterState - Categories, category value and their counts with the current filter applied. Required when
 * @param ontologyTermLabelsById - Set of ontology term labels keyed by term ID, used to determine labels for ontology.
 * checking enabled state of view that is dependent on the state of another category.
 * @returns Ontology view model.
 */
function buildCuratedOntologyCategoryView(
  config: CuratedOntologyCategoryFilterConfig,
  categoryValueByValue: KeyedSelectCategoryValue,
  filterState: FilterState,
  ontologyTermLabelsById: Map<string, string>
): OntologyCategoryView {
  const {
    categoryFilterId,
    isLabelVisible,
    isSearchable,
    isZerosVisible,
    label,
    mask,
  } = config;

  // Build tree view models (e.g. individual tree structures for displaying different ontologies (e.g. human vs mouse
  // vs other for development stage, or just tissues for tissue).
  const treeViews = Object.keys(mask).reduce(
    (accum, ontologyViewKey: string) => {
      const ontologyNodes = mask[ontologyViewKey as ONTOLOGY_VIEW_KEY];
      if (!ontologyNodes) {
        return accum; // Error state - ignore species view.
      }

      // Handle special cases where species is to be excluded.
      if (
        categoryFilterId === CATEGORY_FILTER_ID.DEVELOPMENT_STAGE &&
        !isDevelopmentStageSpeciesVisible(
          filterState,
          ontologyViewKey as ONTOLOGY_VIEW_KEY
        )
      ) {
        return accum;
      }

      // Build view model for each node.
      const childrenViews = ontologyNodes.map((ontologyNode) =>
        buildCuratedOntologyCategoryValueView(
          ontologyNode,
          categoryValueByValue,
          ontologyTermLabelsById
        )
      );

      // Calculate partial selected states.
      childrenViews.forEach((childView) =>
        markOntologySelectedPartialViews(childView)
      );

      // Determine the set of selected values that should be included as tags. Adding this as a flat set to prevent
      // UI code from having to recursively determine selected set when rendering tags.
      const selectedViews = childrenViews.reduce((accum, childView) => {
        listOntologySelectedViews(childView, accum);
        return accum;
      }, new Set<OntologyCategoryTreeNodeView>());

      const viewLabel =
        ONTOLOGY_VIEW_LABEL[
          ontologyViewKey as keyof typeof ONTOLOGY_VIEW_LABEL
        ];

      accum.push({
        children: childrenViews,
        label: isLabelVisible ? viewLabel : undefined,
        selectedViews: [...selectedViews.values()],
      });

      return accum;
    },
    [] as OntologyCategoryTreeView[]
  );

  // Build up the ontology category view model.
  const ontologyView: OntologyCategoryView = {
    isSearchable,
    isZerosVisible,
    key: categoryFilterId,
    label,
    views: treeViews,
  };

  // Check if ontology category is disabled.
  if (isOntologyCategoryViewDisabled(ontologyView)) {
    ontologyView.isDisabled = true;
    ontologyView.tooltip = TOOLTIP_CATEGORY_DISABLED;
  }

  return ontologyView;
}

/**
 * Build view model of node of ontology tree to be displayed as a value in an ontology menu. This is specific to curated
 * ontology category filters such as development stage.
 * @param ontologyNode - Ontology node to build view model for.
 * @param categoryValueByValue - Internal filter model of ontology category.
 * @param ontologyTermLabelsById - Set of ontology term labels keyed by term ID, used to determine labels for ontology.
 * @returns Ontology view model.
 */
function buildCuratedOntologyCategoryValueView(
  ontologyNode: OntologyNode,
  categoryValueByValue: KeyedSelectCategoryValue,
  ontologyTermLabelsById: Map<string, string>
): OntologyCategoryTreeNodeView {
  const { ontology_term_id: categoryValueKey } = ontologyNode;
  const categoryValue = categoryValueByValue.get(categoryValueKey);

  // If there's no corresponding category for this node, create a basic model with 0 count and unselected.
  if (!categoryValue) {
    return {
      count: 0,
      key: categoryValueKey,
      label: ontologyNode.label,
      selected: false,
      selectedPartial: false,
      value: categoryValueKey,
    };
  }

  // Build up base view model.
  const view = {
    count: categoryValue.count,
    key: categoryValueKey,
    label: ontologyTermLabelsById.get(categoryValueKey) ?? categoryValueKey,
    value: categoryValueKey,
  };

  // If this ontology node is a leaf, add its selected value and return.
  if (!ontologyNode.children) {
    return {
      ...view,
      selected: categoryValue.selected,
      selectedPartial: false,
    };
  }

  // Otherwise, build view models for child nodes.
  const children = ontologyNode.children.map((childNode) =>
    buildCuratedOntologyCategoryValueView(
      childNode,
      categoryValueByValue,
      ontologyTermLabelsById
    )
  );

  // Build up view model for this node, including children.
  return {
    ...view,
    children,
    selected: categoryValue.selected,
    selectedPartial: false,
  };
}

// TODO(cc)
interface OntologyPanelCategoryViewBuilder {
  panel: CategoryFilterPanelConfig;
  selectCategoryValues: SelectCategoryValue[];
  selectedValues: CategoryValueId[];
}

/**
 * Build view model of multi-panel category.
 * @param config - Config model of ontology category.
 * @param categoryValueByValue - Internal filter model of single or multiselect category.
 * checking enabled state of view that is dependent on the state of another category.
 * @param multiPanelUIState - Current set of category values that the user has selected in multi-panel cateogry filters.
 * @param ontologyTermLabelsById - Set of ontology term labels keyed by term ID, used to determine labels for ontology
 * @returns Select category view model.
 */
function buildMultiPanelCategoryView(
  config: OntologyMultiPanelFilterConfig,
  categoryValueByValue: KeyedSelectCategoryValue,
  multiPanelUIState: MultiPanelUIState,
  ontologyTermLabelsById: Map<string, string>
): OntologyMultiPanelCategoryView {
  const { categoryFilterId, panels: panelConfigs } = config;
  const categoryFilterUIState = multiPanelUIState.get(categoryFilterId);
  if (!categoryFilterUIState) {
    console.log(
      `Multi-panel category filter state not found for category ${categoryFilterId}`
    );
    return {
      key: categoryFilterId,
      label: config.label,
      panels: [],
      selectedViews: [],
    };
  }

  // Build builders for each panel. TODO(cc) remove now that we have ui state?
  const selectCategoryValues = [...categoryValueByValue.values()];
  const selectedValues = categoryFilterUIState.selected;
  const builders = buildParentPanelBuilders(
    panelConfigs,
    selectCategoryValues,
    selectedValues
  );

  // Build value view models for each panel.
  const ontologyPanelCategoryViews = builders.reduce(
    (accum, builder, index) => {
      // Determine the set of selected values that could possibly restrict the include list for this panel. That is,
      // find the selected values of panels that are parents to this panel.
      const parentBuilders = builders.slice(0, index);
      const visibleSelectCategoryValues = parentBuilders.length
        ? applySelectCategoryValueIncludeList(parentBuilders, builder)
        : builder.selectCategoryValues;

      // Build up view models for each select category value.
      const panelSelectCategoryValueViews: SelectCategoryValueView[] =
        buildSelectCategoryValueViews(
          config,
          visibleSelectCategoryValues,
          ontologyTermLabelsById
        );

      // Apply non-specific labels to values that appear in this panel as well as parent panels. TODO(cc) add discriminating union for this as we only want it for tissue
      const parentPanelView = accum.slice(0, index);
      if (parentPanelView.length) {
        applyNonSpecificLabel(parentPanelView, panelSelectCategoryValueViews);
      }

      // Sort views.
      panelSelectCategoryValueViews.sort(sortCategoryValueViews("label"));

      // Build panel view.
      accum.push({
        label: builder.panel.label,
        views: panelSelectCategoryValueViews,
      });
      return accum;
    },
    [] as OntologyPanelCategoryView[]
  );

  // Determine set of selected values for this multi-panel category
  const allCategoryValueViews = ontologyPanelCategoryViews
    .map((ontologyPanelCategoryView) => ontologyPanelCategoryView.views)
    .flat();
  const selectedViews = buildSelectedViews(
    allCategoryValueViews,
    categoryFilterUIState
  );

  // Build view model of multi-panel category.
  return {
    key: categoryFilterId,
    label: config.label,
    panels: ontologyPanelCategoryViews,
    selectedViews: [...selectedViews],
  };
}

/**
 * TODO(cc) on rename, rename tests, can we run this off just categoryFilterUIState (and just pull the labels from keyed categoryValueViews)
 *
 */
export function buildSelectedViews(
  categoryValueViews: SelectCategoryValueView[],
  categoryFilterUIState: MultiPanelCategoryFilterUIState
): SelectCategoryValueView[] {
  const { selected, selectedPartial, uiNodesByCategoryValueId } =
    categoryFilterUIState;
  // Check if we can add any views to the select set, used to display selected tags.
  return categoryValueViews
    .filter((categoryValueView) => categoryValueView.selected)
    .reduce((accum, categoryValueView) => {
      if (
        isSelectedViewTagVisible(
          categoryValueView.key,
          selected,
          selectedPartial,
          uiNodesByCategoryValueId
        )
      ) {
        accum.push(categoryValueView);
      }
      return accum;
    }, [] as SelectCategoryValueView[]);
}

/**
 * TODO(cc) location, docs, name
 */
function isSelectedViewTagVisible(
  categoryValueId: CategoryValueId,
  selected: CategoryValueId[],
  selectedPartial: CategoryValueId[],
  uiNodesByCategoryValueId: Map<CategoryValueId, MultiPanelUINode>
): boolean {
  // Grab the parents for this selected view.
  const uiParents = getUIParents(categoryValueId, uiNodesByCategoryValueId);

  // If there are no parents, add the selected view.
  if (!uiParents.length) {
    return true;
  }

  // Grab the selected parents.
  const selectedParents = uiParents.filter((uiParent) =>
    selected.includes(uiParent)
  );
  if (!selectedParents.length) {
    // If there are no selected parents, check all ancestors to see if they are included. For example, hema system,
    // blood, non-specific, umbilical cord blood, venous blood.
    return uiParents.every((uiParent) =>
      isSelectedViewTagVisible(
        uiParent,
        selected,
        selectedPartial,
        uiNodesByCategoryValueId
      )
    );
  }

  // If any parent is partially selected, add the selected view.
  return selectedParents.some((uiParent) => selectedPartial.includes(uiParent));
}

/**
 * TODO(cc) docs, location, name
 */
function getUIParents(
  categoryValueId: CategoryValueId,
  uiNodesByCategoryValueId: Map<CategoryValueId, MultiPanelUINode>
): CategoryValueId[] {
  return uiNodesByCategoryValueId.get(categoryValueId)?.uiParents ?? [];
}

/**
 * TODO(cc) docs, location, name
 */
function buildParentPanelBuilders(
  panelConfigs: CategoryFilterPanelConfig[],
  selectCategoryValues: SelectCategoryValue[],
  selectedValues: CategoryValueId[]
): OntologyPanelCategoryViewBuilder[] {
  return panelConfigs.reduce((accum, panelConfig) => {
    const panelSelectCategoryValues =
      panelConfig.sourceKind === "CURATED_CATEGORIES"
        ? maskInferredCuratedCategories(panelConfig.mask, selectCategoryValues)
        : maskAllExact(selectCategoryValues);

    // Determine the set of selected values for this panel.
    const panelSelectedValues = selectedValues.filter(
      (selectedValue: CategoryValueId) =>
        panelSelectCategoryValues.some(
          (selectCategoryValue: SelectCategoryValue) =>
            selectCategoryValue.key === selectedValue
        )
    );

    // Add the builder to the set.
    accum.push({
      panel: panelConfig,
      selectCategoryValues: panelSelectCategoryValues,
      selectedValues: panelSelectedValues,
    });
    return accum;
  }, [] as OntologyPanelCategoryViewBuilder[]);
}

/**
 * Update label of values that also appear in parent panels.
 * @param parentPaneView - Parent panel views.
 * @param panelCategoryValueViews - Views to display in the current panel.
 */
function applyNonSpecificLabel(
  parentPaneView: OntologyPanelCategoryView[],
  panelCategoryValueViews: SelectCategoryValueView[]
) {
  // Get the set of all labels in the parent panels. We are looking to distinguish between label collisions.
  const parentLabels = parentPaneView.reduce(
    (accum: CategoryValueId[], parentPaneView) => {
      parentPaneView.views.forEach((view) => accum.push(view.label));
      return accum;
    },
    [] as CategoryValueId[]
  );

  // TODO(cc) immutable
  panelCategoryValueViews.forEach((panelCategoryValueView) => {
    const { label } = panelCategoryValueView;
    if (parentLabels.includes(label)) {
      panelCategoryValueView.label = `${label}${LABEL_SUFFIX_NON_SPECIFIC}`;
    }
  });
}

/**
 * TODO(cc) docs, location
 * @param parentBuilders - "Intermediate" builder objects of parent panels
 * @param builder - Builder object for the current panel.
 */
function applySelectCategoryValueIncludeList(
  parentBuilders: OntologyPanelCategoryViewBuilder[],
  builder: OntologyPanelCategoryViewBuilder
): SelectCategoryValue[] {
  // If there are no selected values, then all select category values can be included for display.
  const allSelectedValues = parentBuilders
    .map((builder) => builder.selectedValues)
    .flat();
  if (!allSelectedValues.length) {
    return builder.selectCategoryValues;
  }

  const selectedValues = overrideSelectedParents(allSelectedValues);

  // Otherwise, only include values that are either selected or descendants of selected values, or if the value itself
  // is selected.
  const includeOntologyTermIds = selectedValues.reduce(
    (accum: string[], selectedValue: string) => {
      const ontologyTermId = removeOntologyTermIdPrefix(selectedValue);
      accum.push(ontologyTermId);
      accum.push(...(TISSUE_DESCENDANTS[ontologyTermId] ?? []));
      return accum;
    },
    []
  );
  return builder.selectCategoryValues.filter((selectCategoryValue) => {
    const ontologyTermId = removeOntologyTermIdPrefix(selectCategoryValue.key);
    return (
      selectCategoryValue.selected ||
      includeOntologyTermIds.includes(ontologyTermId)
    );
  });
}

/**
 * TODO(cc) docs, location, name
 */
function maskAllExact(
  selectCategoryValues: SelectCategoryValue[]
): SelectCategoryValue[] {
  return selectCategoryValues.filter((selectCategoryValue) => {
    return (
      removeOntologyTermId(selectCategoryValue.key) === OrFilterPrefix.EXPLICIT
    );
  });
}

/**
 * TODO(cc) docs, location, name
 */
function maskInferredCuratedCategories(
  mask: OntologyTermSet,
  selectCategoryValues: SelectCategoryValue[]
): SelectCategoryValue[] {
  const curatedOntologyTermIds = [...listOntologyTreeIds(mask)];
  return selectCategoryValues.filter((selectCategoryValue) => {
    const [prefix, ontologyTermId] = splitOntologyTermIdAndPrefix(
      selectCategoryValue.key
    );
    return (
      prefix === OrFilterPrefix.INFERRED &&
      curatedOntologyTermIds.includes(ontologyTermId)
    );
  });
}

/**
 * Build view model of range category.
 * @param config - Config model of ontology category.
 * @param rangeCategory - Internal filter model of range category.
 * @returns Range view model.
 */
function buildRangeCategoryView(
  config: CategoryFilterConfig,
  rangeCategory: RangeCategory
): RangeCategoryView {
  const { categoryFilterId } = config;

  // Build view model of range category.
  const rangeView: RangeCategoryView = {
    key: categoryFilterId,
    label: config.label,
    max: rangeCategory.max,
    min: rangeCategory.min,
    selectedMax: rangeCategory.selectedMax,
    selectedMin: rangeCategory.selectedMin,
  };

  // Determine if range view is disabled.
  if (isRangeCategoryDisabled(rangeView)) {
    rangeView.isDisabled = true;
    rangeView.tooltip = TOOLTIP_CATEGORY_DISABLED;
  }

  return rangeView;
}

/**
 * Build select category value view models of the given select category values.
 * @param config - Config model of category to build category value views for.
 * @param selectCategoryValues - Values to build view models from.
 * @param ontologyTermLabelsById - Set of ontology term labels keyed by term ID, used to determine labels for ontology
 */
function buildSelectCategoryValueViews(
  config: CategoryFilterConfig,
  selectCategoryValues: SelectCategoryValue[],
  ontologyTermLabelsById: Map<string, string>
): SelectCategoryValueView[] {
  return selectCategoryValues.map(
    ({ count, key, selected, selectedPartial }: SelectCategoryValue) => {
      return {
        count,
        key,
        label: buildCategoryValueLabel(config, key, ontologyTermLabelsById),
        selected,
        selectedPartial,
        value: key,
      };
    }
  );
}

/**
 * Build view model of single or multiselect category.
 * @param config - Config model of ontology category.
 * @param categoryValueByValue - Internal filter model of single or multiselect category.
 * @param filterState - Categories, category value and their counts with the current filter applied. Required when
 * @param ontologyTermLabelsById - Set of ontology term labels keyed by term ID, used to determine labels for ontology
 * checking enabled state of view that is dependent on the state of another category.
 * @returns Select category view model.
 */
function buildSelectCategoryView(
  config: CategoryFilterConfig,
  categoryValueByValue: KeyedSelectCategoryValue,
  filterState: FilterState,
  ontologyTermLabelsById: Map<string, string>
): SelectCategoryView {
  const { categoryFilterId } = config;

  // Grab the config for this category.
  const { pinnedCategoryValues, tooltip } =
    CATEGORY_FILTER_CONFIGS_BY_ID[categoryFilterId];

  const allCategoryValueViews = buildSelectCategoryValueViews(
    config,
    [...categoryValueByValue.values()],
    ontologyTermLabelsById
  ).sort(sortCategoryValueViews("key"));

  // Split values into pinned and unpinned.
  const [pinnedValues, unpinnedValues] = partitionSelectCategoryValueViews(
    allCategoryValueViews,
    pinnedCategoryValues
  );

  // Build view model of select category.
  const selectView: SelectCategoryView = {
    key: categoryFilterId,
    label: config.label,
    pinnedValues,
    unpinnedValues,
    values: allCategoryValueViews,
  };

  // Handle special cases where select category may be disabled.
  if (
    categoryFilterId === CATEGORY_FILTER_ID.ETHNICITY &&
    !isEthnicityViewEnabled(filterState)
  ) {
    selectView.isDisabled = true;
    selectView.tooltip = tooltip;
  }
  // Otherwise, check generic case where category is disabled due to no values meeting current filter.
  else if (isSelectCategoryDisabled(selectView)) {
    selectView.isDisabled = true;
    selectView.tooltip = TOOLTIP_CATEGORY_DISABLED;
  }

  return selectView;
}

/**
 * Find and return the selected values for the given category.
 * @param categoryFilterId - ID of category to find selected filters of.
 * @param filters - Current set of selected category values (values) keyed by category (id).
 * @returns Array of filters
 */
function getCategoryFilter<T extends Categories>(
  categoryFilterId: CATEGORY_FILTER_ID,
  filters: Filters<T>
): CategoryFilter | undefined {
  return filters.find((filter) => filter.id === categoryFilterId);
}

/**
 * Reevaluate selected state of parent. It's possible all children of this parent are now selected; if so, add parent
 * to set of selected values and then reevaluate the parent's parent selected state.
 * @param ontologyRootNodes - Top-level nodes in ontology tree.
 * @param ontologyNode - Node to reevaluate selected state of.
 * @param selectedCategoryValues - The current set of selected values.
 * @param categoryKeyValues - Original, full set of values for this category.
 */
function handleOntologyChildAdded(
  ontologyRootNodes: OntologyNode[],
  ontologyNode: OntologyNode,
  selectedCategoryValues: Set<CategoryValueId>,
  categoryKeyValues: Set<CategoryValueId>
) {
  // Check if all children of the parent node are selected, only including children values that are present in the
  // original result set.
  const isEveryChildSelected = isEveryChildNodeSelected(
    ontologyNode,
    selectedCategoryValues,
    categoryKeyValues
  );

  // Add parent if all children are selected.
  if (isEveryChildSelected) {
    selectedCategoryValues.add(ontologyNode.ontology_term_id);

    // Node's parent, if any, must now be reevaluated.
    const parentNode = findOntologyParentNode(ontologyRootNodes, ontologyNode);
    if (parentNode) {
      handleOntologyChildAdded(
        ontologyRootNodes,
        parentNode,
        selectedCategoryValues,
        categoryKeyValues
      );
    }
  }
}

/**
 * Reevaluate selected state of node. If a child has been removed, then a parent can no longer be considered selected;
 * remove node from selected set and reevaluate the parent's parent selected state.
 * @param ontologyRootNodes - Top-level nodes in ontology tree.
 * @param ontologyNode - Node to reevaluate selected state of.
 * @param selectedCategoryValues - The current set of selected values.
 */
function handleOntologyChildRemoved(
  ontologyRootNodes: OntologyNode[],
  ontologyNode: OntologyNode,
  selectedCategoryValues: Set<CategoryValueId>
) {
  const { ontology_term_id: parentId } = ontologyNode;
  if (selectedCategoryValues.has(parentId)) {
    // Parent is currently selected, remove it from the selected set.
    selectedCategoryValues.delete(parentId);

    // Node's parent, if any, must now be reevaluated.
    const parentNode = findOntologyParentNode(ontologyRootNodes, ontologyNode);
    if (parentNode) {
      handleOntologyChildRemoved(
        ontologyRootNodes,
        parentNode,
        selectedCategoryValues
      );
    }
  }
}

/**
 * Determine the rows that have values matching the given filters. Row must have at least one selected value across each
 * category. Mimics react-query's includeSome functionality.
 * @param rowValue - Value to filter row by.
 * @param filter - Selected filter to apply to row.
 * @returns Filtered array of rows.
 */
function includesSome(
  rowValue: string | string[],
  filter: CategoryFilter
): boolean {
  // "or" across category values (that is, inside a category).
  return (
    rowValue &&
    rowValue.length &&
    filter.value.some((val: string) => rowValue.includes(val)) // Handles string or array values
  );
}

/**
 * Check if the given category filter's match (filter) type is "between".
 * @param categoryFilterId - ID of category to check type of.
 * @returns True if the given category's type is "between".
 */
function isMatchKindBetween(categoryFilterId: CATEGORY_FILTER_ID): boolean {
  const config = CATEGORY_FILTER_CONFIGS_BY_ID[categoryFilterId];
  return config.matchKind === "BETWEEN";
}

/**
 * Determine if the given ethnicity is considered unspecified (that is, na or unknown).
 * @param categoryValueKey - Ethnicity value to check if it's specified.
 * @returns True if ethnicity is either na or unknown.
 */

function isEthnicitySpecified(categoryValueKey: CategoryValueId) {
  const label =
    ETHNICITY_UNSPECIFIED_LABEL[
      categoryValueKey as keyof typeof ETHNICITY_UNSPECIFIED_LABEL
    ];
  return !label;
}

/*
 * Check if all children of the parent node are selected, only including children values that are present in the
 * original result set.
 * @param ontologyNode - Node to check selected state of children.
 * @param selectedCategoryValues - The current set of selected values.
 * @param categoryKeyValues - Original, full set of values for this category.
 */
function isEveryChildNodeSelected(
  ontologyNode: OntologyNode,
  selectedCategoryValues: Set<CategoryValueId>,
  categoryKeyValues: Set<CategoryValueId>
): boolean {
  // Check if all children of the parent node are selected, only including children values that are present in the
  // original result set.
  const childrenIds =
    ontologyNode.children
      ?.map((child) => child.ontology_term_id)
      .filter((childId) => categoryKeyValues.has(childId)) ?? [];
  return Boolean(
    childrenIds.length &&
      childrenIds.every((childId) => selectedCategoryValues.has(childId))
  );
}

/**
 * The ethnicity filter is only view enabled if:
 * 1. Homo sapiens is available as an option in the organism filter.
 * 2. The organism filter has selected values that includes Homo sapiens.
 * @param filterState - Categories, category value and their counts with the current filter applied. Required to
 * determine if ethnicity category should be enabled.
 * @returns True if ethnicity is either na or unknown.
 */
function isEthnicityViewEnabled(filterState: FilterState) {
  // Check to see if there are Homo sapiens values in the result set.
  const organismCategoryValues = filterState[
    CATEGORY_FILTER_ID.ORGANISM
  ] as KeyedSelectCategoryValue;
  const count = organismCategoryValues.get(ORGANISM.HOMO_SAPIENS)?.count ?? 0;
  if (count === 0) {
    return false;
  }

  // Check to see if there are any selected values for organism and if so, Homo sapiens must be one of them.
  const selectedOrganisms = [...organismCategoryValues.values()]
    .filter((selectCategoryValue) => selectCategoryValue.selected)
    .map((selectCategoryValue) => selectCategoryValue.key);
  return (
    selectedOrganisms.length === 0 ||
    selectedOrganisms.includes(ORGANISM.HOMO_SAPIENS)
  );
}

/**
 * Determine if given filters are identical.
 * @param filters0 - First filter to compare.
 * @param filters1 - Second filter to compare.
 */
function isFilterEqual<T extends Categories>(
  filters0: Filters<T>,
  filters1: Filters<T>
): boolean {
  return (
    filters0.length === filters1.length &&
    filters0.every((val) => filters1.includes(val))
  );
}

/**
 * Determine if the given selected value is a range and a category value ID.
 * @param selectedValue - Selected filter value, either a category value ID (e.g. "normal"), or a
 * range (e.g. [1,100]).
 * @returns True if given selected value is a range.
 */
function isSelectedValueRange(
  selectedValue: CategoryValueId | Range
): selectedValue is Range {
  return Array.isArray(selectedValue);
}

/**
 * Determine if the given category value is a select category value (and not a range category value).
 * @param categorySetValue - Range or category set value.
 * @returns True if category set value is a set of category value keys.
 */
function isCategoryValueIdSet(
  categorySetValue: CategorySetValue
): categorySetValue is Set<CategoryValueId> {
  return categorySetValue instanceof Set;
}

/**
 * Determine if the given category config is for a curated ontology category filter.
 * @param config - Config model of category, either an ontology category config or a base category config.
 * @returns True if category config is for a curated ontology category.
 */
function isCuratedOntologyCategoryFilterConfig(
  config: CategoryFilterConfig
): config is CuratedOntologyCategoryFilterConfig {
  return config.valueSourceKind == "CURATED";
}

/**
 * Determine if the given category config is a multi-panel category.
 * @param config - Config model of category, either an ontology category config or a base category config.
 * @returns True if category config is for a multi-panel category.
 */
function isMultiPanelCategoryFilterConfig(
  config: CategoryFilterConfig
): config is OntologyMultiPanelFilterConfig {
  return config.viewKind === "MULTI_PANEL";
}

/**
 * Returns true if ontology category is disabled, that is, there are no species ontology trees or species
 * ontology trees with values that have a count.
 * @param categoryView - Ontology category view to check enabled/disabled state of.
 * @returns True when there are no species or the count of children for each species is 0.
 */
function isOntologyCategoryViewDisabled(
  categoryView: OntologyCategoryView
): boolean {
  const { views } = categoryView;
  if (!views || views.length === 0) {
    return true;
  }
  return !views.some((s) => s.children.some((child) => child.count > 0));
}

/**
 * Returns true if range category is disabled, that is, range min and max are both 0 or both equal.
 * @param categoryView - Range category view to check enabled/disabled state of.
 * @returns true when range min and max are both 0 or both equal.
 */
function isRangeCategoryDisabled(categoryView: RangeCategoryView): boolean {
  const { max, min } = categoryView;
  return (min === 0 && max === 0) || min === max;
}

/**
 * Returns true if select category is disabled, that is, the category is disabled or all values have a count of 0.
 * @param categoryView - Select category view to check enabled/disabled state of.
 * @returns true when the category is disabled or all category values have a count of 0.
 */
function isSelectCategoryDisabled(categoryView: SelectCategoryView): boolean {
  const { isDisabled, values } = categoryView;
  return isDisabled || values.every((value) => value.count === 0);
}

/**
 * Determine if the given category value is a select (including ontology) category value (and not a range category value).
 * @param categoryValue - Range or select (including ontology) category value.
 * @returns True if category value is a select category value.
 */
function isSelectCategoryValue(
  categoryValue: RangeCategory | KeyedSelectCategoryValue
): categoryValue is KeyedSelectCategoryValue {
  return categoryValue instanceof Map;
}

/**
 * Development stage species is only visible if:
 * 1. There are no selected organisms or,
 * 2. The given species is selected.
 * @param filterState - Categories, category value and their counts with the current filter applied. Required to
 * determine if development stfage species should be visible.
 * @param speciesKey - The species to check if a corresponding organism has been selected for.
 * @returns True if given species is to be displayed.
 */
function isDevelopmentStageSpeciesVisible(
  filterState: FilterState,
  speciesKey: ONTOLOGY_VIEW_KEY
) {
  // Find the current selected values for organism.
  const organismCategoryValues = filterState[
    CATEGORY_FILTER_ID.ORGANISM
  ] as KeyedSelectCategoryValue;
  const selectedOrganisms = [...organismCategoryValues.values()]
    .filter((selectCategoryValue) => selectCategoryValue.selected)
    .map((selectCategoryValue) => selectCategoryValue.key);

  // If no organisms are selected, all species can be displayed.
  if (selectedOrganisms.length === 0) {
    return true;
  }

  // Otherwise this species is only visible if it's selected.
  if (speciesKey === ONTOLOGY_VIEW_KEY.HsapDv) {
    return selectedOrganisms.includes(ORGANISM.HOMO_SAPIENS);
  }
  if (speciesKey === ONTOLOGY_VIEW_KEY.MmusDv) {
    return selectedOrganisms.includes(ORGANISM.MUS_MUSCULUS);
  }
  // Check the "other" case where any species other than human and mouse must be selected.
  return (
    selectedOrganisms.filter(
      (organism) =>
        organism !== ORGANISM.HOMO_SAPIENS && organism !== ORGANISM.MUS_MUSCULUS
    ).length > 0
  );
}

/**
 * Build up model of ontology term labels keyed by ontology term ID for the given ontology nodes. This map is used
 * to generate display labels for category filters that are backed by ontology term IDs.
 * @param labelsById - Map of ontology term labels keyed by ontology term ID.
 * @param ontologyNodes - Array of ontology nodes to key labels for.
 */
function keyCuratedOntologyNodeTermLabelsById(
  labelsById: Map<string, string>,
  ontologyNodes: OntologyNode[]
) {
  ontologyNodes.forEach((ontologyNode) => {
    if (!labelsById.has(ontologyNode.label)) {
      labelsById.set(ontologyNode.ontology_term_id, ontologyNode.label);
      if (ontologyNode.children) {
        keyCuratedOntologyNodeTermLabelsById(labelsById, ontologyNode.children);
      }
    }
  });
}

/**
 * Build up model of ontology term labels keyed by ontology term ID for the given ontology term set. This map is
 * used to generate display labels for category filters that are backed by ontology term IDs.
 * @param labelsById - Map of ontology term labels keyed by ontology term ID.
 * @param mask - Set of ontology terms to key labels for.
 */
function keyCuratedOntologyTermLabelsById(
  labelsById: Map<string, string>,
  mask: OntologyTermSet
) {
  Object.keys(mask).forEach((ontologyViewKey: string) => {
    const ontologyNodes = mask[ontologyViewKey as ONTOLOGY_VIEW_KEY];
    if (!ontologyNodes) {
      return labelsById; // Error state - ignore species view.
    }
    keyCuratedOntologyNodeTermLabelsById(labelsById, ontologyNodes);
  });
}

/**
 * Build up model of ontology term labels keyed by ontology term ID for the given set of rows and row value. This map is
 * used to generate display labels for category filters that are backed by ontology term IDs.
 * @param labelsById - Map of ontology term labels keyed by ontology term ID.
 * @param rows - Original result set to key ontology term labels for.
 * @param filterOnKey - Row value key to key ontology term labels for.
 */
function keyRowOntologyTermLabelsById<T extends Categories>(
  labelsById: Map<string, string>,
  rows: Row<T>[],
  filterOnKey: FilterKey
) {
  rows.forEach((originalRow: Row<T>) => {
    // eslint-disable-next-line @typescript-eslint/ban-ts-comment --- TODO(cc) revisit - different between FilterKey (any value from dataset or collection) vs T extends Categories, also see if type assertion can be resolved
    // @ts-ignore -- as above
    (originalRow.original[filterOnKey] as Ontology[]).forEach(
      ({ label, ontology_term_id }) => {
        if (!labelsById.has(ontology_term_id)) {
          labelsById.set(ontology_term_id, label);
        }
      }
    );
  });
}

/**
 * Build map of ontology term labels keyed by ontology term IDs. This map is used to determine the display values for
 * ontology term-backed fields (e.g. tissue).
 * @param originalRows - Original result set before filtering.
 * @returns Map of ontology term labels keyed by ID.
 * TODO(cc) revisit this - need to do both top-level configs and panel configs. can we make this configurable?
 */
function keyOntologyTermLabelsById<T extends Categories>(
  originalRows: Row<T>[]
): Map<string, string> {
  const labelsById = new Map<string, string>();

  // Key curated term sets.
  [
    DEVELOPMENT_STAGE_ONTOLOGY_TERM_SET,
    TISSUE_SYSTEM_ONTOLOGY_TERM_SET,
    TISSUE_ORGAN_ONTOLOGY_TERM_SET,
  ].forEach((termSet: OntologyTermSet) => {
    keyCuratedOntologyTermLabelsById(labelsById, termSet);
  });

  // Key ontology terms from row values.
  ["cell_type", "tissue"].forEach((filterOnKey) => {
    keyRowOntologyTermLabelsById(
      labelsById,
      originalRows,
      filterOnKey as FilterKey
    );
  });

  return labelsById;
}

/**
 * Update partially selected state of views. If parent isn't selected but some of its children are selected or partially
 * selected, mark parent as partially selected
 * @param view - View model of ontology category value.
 */
function markOntologySelectedPartialViews(view: OntologyCategoryTreeNodeView) {
  // Mark any children as partially selected.
  view.children?.forEach((childView) => {
    markOntologySelectedPartialViews(childView);
  });

  // If this view isn't selected, check if any of it's children are either selected or partially selected.
  if (view.selected) {
    return;
  }
  view.selectedPartial =
    view.children?.some(
      (childView) => childView.selected || childView.selectedPartial
    ) ?? false;
}

/**
 * Determine the set of selected values (views) that form the basis of the set of selected tags. Adding this as a flat
 * set to prevent UI code from having to recursively determine selected set.
 * @param view - View model of ontology category value.
 * @param selectedSet - Selected set of cateogry value keys.
 */
function listOntologySelectedViews(
  view: OntologyCategoryTreeNodeView,
  selectedSet: Set<OntologyCategoryTreeNodeView>
) {
  // If the view is selected it can be included in the selected set of tags.
  if (view.selected) {
    selectedSet.add(view);
    return;
  }

  // If the view is partially selected, check its children to see which ones can be included in the set of selected tags.
  if (view.selectedPartial && view.children) {
    view.children.forEach((childView) =>
      listOntologySelectedViews(childView, selectedSet)
    );
  }

  return selectedSet;
}

/**
 * Split select category values into arrays of pinned and non-pinned values.
 * @param categoryValues - Category value view models for a given category.
 * @param pinnedCategoryValues - Category value keys to be pinned for a given category.
 * @returns Tuple containing an array of pinned values and an array of non-pinned values.
 */
function partitionSelectCategoryValueViews(
  categoryValues: SelectCategoryValueView[],
  pinnedCategoryValues?: CategoryValueId[]
): [SelectCategoryValueView[], SelectCategoryValueView[]] {
  // Handle case where category has no pinned values.
  if (!pinnedCategoryValues) {
    return [[], categoryValues];
  }

  // Otherwise, split category values into pinned and non-pinned arrays.
  const partitionedValues: [
    SelectCategoryValueView[],
    SelectCategoryValueView[]
  ] = [[], []];
  return categoryValues.reduce((accum, categoryValue) => {
    const [pinned, nonPinned] = accum;
    if (pinnedCategoryValues.includes(categoryValue.key)) {
      pinned.push(categoryValue);
    } else {
      nonPinned.push(categoryValue);
    }
    return accum;
  }, partitionedValues);
}

/**
 * Remove all descendents from the set of selected values.
 * @param ontologyNodes - Nodes to remove from selected set of values.
 * @param selectedCategoryValues - The current set of selected values.
 */
function removeOntologyDescendents(
  ontologyNodes: OntologyNode[],
  selectedCategoryValues: Set<CategoryValueId>
) {
  // For each node, remove self and children.
  ontologyNodes.forEach((ontologyNode) => {
    selectedCategoryValues.delete(ontologyNode.ontology_term_id);
    if (ontologyNode.children) {
      removeOntologyDescendents(ontologyNode.children, selectedCategoryValues);
    }
  });
}

/**
 * Handle select of multi-panel value: build and set next set of filters for this category. Track selected select value.
 * @param config - Configuration model of selected category.
 * @param categoryValueId - The selected category value.
 * @param selectedValue - Selected category value ID to use as selected value.
 * @param source - TODO(cc)
 * @param setFilter - Function to update set of selected values for a category.
 * @param multiPanelUIState - Current set of category values that the user has selected on the UI for all multi-
 * panel category filters.
 * @param setMultiPanelUIState - React state mutator for setting multi-panel UI state.
 */
// eslint-disable-next-line sonarjs/cognitive-complexity -- TODO(cc) revisit
function onFilterMultiPanelCategory(
  config: OntologyMultiPanelFilterConfig,
  categoryValueId: CategoryValueId,
  selectedValue: CategoryValueId,
  source: ON_FILTER_SOURCE,
  setFilter: SetFilterFn,
  multiPanelUIState: MultiPanelUIState,
  setMultiPanelUIState: Dispatch<SetStateAction<MultiPanelUIState>>
) {
  const { categoryFilterId } = config;

  // Grab the selected values for this category filter.
  const multiPanelFilters = buildMultiPanelFilters(multiPanelUIState);

  // Track selected category and value. TODO(cc) revisit categoryFilters here
  trackSelectCategoryValueSelected(config, categoryValueId, multiPanelFilters);

  // Determine selected set of values for this category filter based on the current selected values for this multi-panel
  // category filter; toggle current selected values.
  const allSelectedCategoryFilters = buildNextSelectCategoryFilters(
    config,
    selectedValue,
    multiPanelFilters
  );

  // If the source of the filter action was a selected tag, we know a remove action has occurred. Clear out any
  // children of the removed tag unless the children themselves have another that is partially/selected.
  // TODO(cc) move to function with buildNextSelectCategoryFilters above, have a single return value of allSelectedCategoryFilters then remove ref's to this set
  const filteredSelectedCategoryFilters = new Set<CategoryValueId>(
    allSelectedCategoryFilters
  );
  if (source === ON_FILTER_SOURCE.TAG) {
    // Grab all selected descendants of the removed value.
    const selectedDescendants = allSelectedCategoryFilters.filter(
      (selectedCategoryValueId) =>
        isDescendant(
          selectedCategoryValueId,
          categoryValueId,
          TISSUE_DESCENDANTS
        )
    );

    // Remove the selected descendant unless the descendant has another parent that is partially/selected.
    // TODO(cc) handle undefined's here and optional chaining here
    const categoryFilterUIState =
      multiPanelUIState.get(categoryFilterId) ??
      ({} as MultiPanelCategoryFilterUIState);
    const { selected, selectedPartial, uiNodesByCategoryValueId } =
      categoryFilterUIState;

    selectedDescendants.forEach((selectedDescendant) => {
      const uiParents =
        uiNodesByCategoryValueId.get(selectedDescendant)?.uiParents; // TODO(cc) chaining

      // If there are no parents, remove the descendant from the selected list.
      if (!uiParents || !uiParents.length) {
        filteredSelectedCategoryFilters.delete(selectedDescendant);
      } else {
        // Don't delete selected descendant if any parent is partially/selected
        const isAnyParentSelected = uiParents
          .filter(
            (uiParent) =>
              uiParent !== categoryValueId &&
              !selectedDescendants.includes(uiParent)
          )
          .some(
            (uiParent) =>
              selected.includes(uiParent) || selectedPartial.includes(uiParent)
          );
        if (!isAnyParentSelected) {
          filteredSelectedCategoryFilters.delete(selectedDescendant);
        }
      }
    });
  }

  // Determine the selected values to pass to react-table to filter the rows by applying restrictions across selected
  // filter values. For example, even if "digestive system" and "tongue" are both selected in the UI, we only want to
  // pass the more restrictive "tongue" to react-table as we only want to show rows that match "tongue".
  const overriddenSelectedCategoryFilters = overrideSelectedParents([
    ...filteredSelectedCategoryFilters,
  ]);

  // Determine the selected values which are now partially selected, if any. For example, if "digestive system" and
  // "tongue" are both selected in the UI, internally we record "digestive system" as partially selected and "tongue"
  // as selected.
  const multiPanelSelectedValues = multiPanelUIState.get(categoryFilterId); // TODO(cc) rename
  if (!multiPanelSelectedValues) {
    console.log(
      `Selected values not found for category filter ID ${categoryFilterId}`
    );
    return;
  }
  const { uiNodesByCategoryValueId } = multiPanelSelectedValues;
  const selectedPartial = listPartiallySelectedCategoryValueIds(
    [...filteredSelectedCategoryFilters],
    overriddenSelectedCategoryFilters,
    uiNodesByCategoryValueId
  );

  // Update internal multi-panel filter state with the updated, toggled set of selected filters for this category as
  // well as the partially (overridden) selected values.
  multiPanelUIState.set(categoryFilterId, {
    selected: [...filteredSelectedCategoryFilters],
    selectedPartial,
    uiNodesByCategoryValueId,
  });
  setMultiPanelUIState(multiPanelUIState);

  // Trigger filter of rows in react-table.
  setFilter(categoryFilterId, overriddenSelectedCategoryFilters);
}

/**
 * We need to know are you blocked by some of your children being selected (partially selected), or by all of your
 * children being selected (not partially selected).
 * TODO(cc) location, name, rename tests
 */
export function listPartiallySelectedCategoryValueIds(
  selectedCategoryValueIds: CategoryValueId[],
  overriddenSelectedCategoryValueIds: CategoryValueId[],
  uiNodesByCategoryValueId: Map<CategoryValueId, MultiPanelUINode>
): CategoryValueId[] {
  const selectedPartial: CategoryValueId[] = [];

  selectedCategoryValueIds
    // Ignore values that are in the override list; these values are still selected and don't need to be identified
    // as partially selected.
    .filter(
      (selectedCategoryFilter) =>
        !overriddenSelectedCategoryValueIds.includes(selectedCategoryFilter)
    )
    // Otherwise, value has been blocked by a more "precise" value. If all children of the blocked value are selected,
    // leave as is. If only some children of the blocked value are selected, add selected value to partial list.
    .forEach((selectedCategoryFilter: CategoryValueId) => {
      const uiNode = uiNodesByCategoryValueId.get(selectedCategoryFilter);
      if (!uiNode) {
        // TODO(cc) is this even possible?
        return;
      }

      const isEveryChildSelected = uiNode.uiChildren.every((uiChild) => {
        // Confirm child is selected and not overridden.
        const isChildSelected =
          selectedCategoryValueIds.includes(uiChild) &&
          overriddenSelectedCategoryValueIds.includes(uiChild);
        if (isChildSelected) {
          return true;
        }

        // If the child is not explicitly selected and it's an exact value, exit here.
        if (isExplicitOntologyTermId(uiChild)) {
          return false;
        }

        // Otherwise, if child is not explicitly selected and it's an inferred value, check if it's inferred selected.
        // TODO(cc) can we move this to a format similar to buildSelectedViews (ie recursive?)
        const uiGrandchildren =
          uiNodesByCategoryValueId.get(uiChild)?.uiChildren;
        if (!uiGrandchildren) {
          return false; // If there are no grandchildren then the value is not inferred selected.
        }
        return uiGrandchildren.every((uiGrandchild) => {
          return selectedCategoryValueIds.includes(uiGrandchild);
        });
      });

      if (!isEveryChildSelected) {
        selectedPartial.push(selectedCategoryFilter);
      }
    });

  return selectedPartial;
}

/**
 * Build up react-table filters model from the given multi-panel UI state.
 * @param multiPanelUIState - ID of category filter to filters.
 * @returns Set of selected values in the format expected by react-table.
 */
function buildMultiPanelFilters<T extends Categories>(
  multiPanelUIState: MultiPanelUIState
): Filters<T> {
  return [...multiPanelUIState.keys()].reduce(
    (accum: CategoryFilter[], categoryFilterId: CATEGORY_FILTER_ID) => {
      // Don't add category to filters if it has no selected values.
      const selectedCategoryValues =
        multiPanelUIState.get(categoryFilterId)?.selected;
      if (!selectedCategoryValues || !selectedCategoryValues.length) {
        return accum;
      }

      // Otherwise, this category has selected values: add them!
      accum.push({
        id: categoryFilterId,
        value: selectedCategoryValues,
      });
      return accum;
    },
    [] as CategoryFilter[]
  );
}

/**
 * TODO(cc) location, rename, can we remove this and use the partial state instead
 */
function overrideSelectedParents(
  selectedValues: CategoryValueId[]
): CategoryValueId[] {
  const selectedOntologyTermIds = selectedValues.map((selectedValue) =>
    removeOntologyTermIdPrefix(selectedValue)
  );
  return selectedValues.reduce((accum, selectedValue: CategoryValueId) => {
    const [selectedPrefix, selectedOntologyId] =
      splitOntologyTermIdAndPrefix(selectedValue);

    // If the selected value is an explicit value, always include it in the selected set of values.
    if (selectedPrefix === OrFilterPrefix.EXPLICIT) {
      accum.push(selectedValue);
      return accum;
    }

    // Handle self case. Check if the corresponding explicit value is selected. For example, "eye" organ (which is
    // inferred) is a parent of "eye" tissue (which is explicit). In this case, the explicit "eye" tissue takes
    // precedence over the "eye" organ.
    const explicitSelectedValue =
      buildExplicitOntologyTermId(selectedOntologyId);
    if (selectedValues.includes(explicitSelectedValue)) {
      return accum;
    }

    // Otherwise, if any descendant of the selected value is selected, do not include the selected value in the set of
    // selected values.
    const descendants = TISSUE_DESCENDANTS[selectedOntologyId] ?? [];

    // Check if any descendant is selected.
    const isAnyDescendantSelected = descendants.some((descendant) =>
      selectedOntologyTermIds.includes(descendant)
    );
    if (!isAnyDescendantSelected) {
      accum.push(selectedValue);
    }

    return accum;
  }, [] as CategoryValueId[]);
}

/**
 * Handle select of ontology value: build and set next set of filters for this category. Track selected ontology value.
 * @param config - Configuration model of selected category.
 * @param categoryValueKey - The selected category value.
 * @param selectedValue - Selected category value ID to use as selected value.
 * @param filters - Current set of selected category values (values) or ranges keyed by category (id).
 * @param setFilter - Function to update set of selected values for a category.
 * @param categorySet - Original, unfiltered sets of category values keyed by their category.
 */
function onFilterCuratedOntologyCategory<T extends Categories>(
  config: CuratedOntologyCategoryFilterConfig,
  categoryValueKey: CategoryValueId,
  selectedValue: CategoryValueId,
  filters: Filters<T>,
  setFilter: SetFilterFn,
  categorySet: CategorySet
) {
  const { categoryFilterId, mask } = config;

  // Track selected category and value.
  trackCuratedOntologyCategoryValueSelected(config, categoryValueKey, filters);

  // Build and set next set of filters for this category.
  const nextCategoryFilters = buildNextOntologyCategoryFilters(
    categoryFilterId,
    selectedValue,
    filters,
    categorySet[categoryFilterId] as Set<CategoryValueId>,
    mask
  );
  setFilter(categoryFilterId, nextCategoryFilters);
}

/**
 * Handle select of range min/max value: set next set of filters for this category. Track updated range.
 * @param config - Configuration model of selected category.
 * @param selectedValue - Selected category value key (e.g. [1, 100]).
 * @param setFilter - Function to update set of selected values for a category.
 */
function onFilterRangeCategory(
  config: CategoryFilterConfig,
  selectedValue: Range,
  setFilter: SetFilterFn
) {
  const { analyticsEvent, categoryFilterId } = config;

  // Track select of new range mim/max, ignoring any clear of selected range. Only track if event is specified on
  // configuration model
  if (analyticsEvent && selectedValue.length > 0) {
    const [min, max] = selectedValue;
    track(analyticsEvent, {
      max,
      min,
    });
  }

  // Update filters for this range category. Convert range min/max object to tuple for react-table.
  setFilter(categoryFilterId, selectedValue);
}

/**
 * Handle select of select value: build and set next set of filters for this category. Track selected select value.
 * @param config - Configuration model of selected category.
 * @param categoryValueKey - The selected category value.
 * @param selectedValue - Category value ID to use as selected value.
 * @param filters - Current set of selected category values (values) or ranges keyed by category (id).
 * @param setFilter - Function to update set of selected values for a category.
 */
function onFilterSelectCategory<T extends Categories>(
  config: CategoryFilterConfig,
  categoryValueKey: CategoryValueId,
  selectedValue: CategoryValueId,
  filters: Filters<T>,
  setFilter: SetFilterFn
) {
  const { categoryFilterId } = config;

  // Track selected category and value.
  trackSelectCategoryValueSelected(config, categoryValueKey, filters);

  // Build and set next set of filters for this category.
  const nextFilters = buildNextSelectCategoryFilters(
    config,
    selectedValue,
    filters
  );

  setFilter(categoryFilterId, nextFilters);
}

/**
 * Update selected state of categories to match the current set of selected filters.
 * @param nextFilterState - Filter state being built on select of filter.
 * @param filters - Current set of selected category values (values) or ranges keyed by category (id).
 * @param multiPanelUIState - Current set of category values that the user has selected on the UI for all multi-
 * panel category filters.
 */
function setSelectedStates<T extends Categories>(
  nextFilterState: FilterState,
  filters: Filters<T>,
  multiPanelUIState: MultiPanelUIState
) {
  Object.keys(nextFilterState).forEach((strCategoryFilterId: string) => {
    const categoryFilterId = strCategoryFilterId as CATEGORY_FILTER_ID;

    // Grab the filter state for this category.
    const categoryFilterState = nextFilterState[categoryFilterId];

    // Grab the filters for this category, used to determine the selected values.
    const config = CATEGORY_FILTER_CONFIGS_BY_ID[categoryFilterId];
    const isMultiPanel = isMultiPanelCategoryFilterConfig(config);
    const categoryFilter = isMultiPanel
      ? getCategoryFilter(
          categoryFilterId,
          buildMultiPanelFilters(multiPanelUIState)
        )
      : getCategoryFilter(categoryFilterId, filters);

    // Grab the partially selected values if this is a multi-panel category.
    const selectedPartials = isMultiPanel
      ? multiPanelUIState.get(categoryFilterId)?.selectedPartial ?? []
      : [];

    if (
      !categoryFilter ||
      (!categoryFilter.value && !selectedPartials.length)
    ) {
      // There are no selected or partially selected values for this category
      return;
    }

    // Handle single and multiselect categories, or ontology categories.
    if (isSelectCategoryValue(categoryFilterState)) {
      // Create sets for easy lookup of category values.
      const selectedCategoryValues = new Set(categoryFilter.value);
      const selectedPartialCategoryValues = new Set(selectedPartials);

      // Check each category value in this category to see if it's selected.
      for (const [
        categoryValueKey,
        categoryValue,
      ] of categoryFilterState.entries()) {
        const selectedPartial =
          selectedPartialCategoryValues.has(categoryValueKey);
        categoryValue.selected =
          selectedCategoryValues.has(categoryValueKey) && !selectedPartial;
        categoryValue.selectedPartial = selectedPartial;
      }
      return;
    }

    // Handle range categories.
    const [selectedMin, selectedMax] = categoryFilter.value;
    categoryFilterState.selectedMin = selectedMin;
    categoryFilterState.selectedMax = selectedMax;
  });
}

/**
 * Sort category value views by the given key, ascending. Sort can be by either the key (required for values such
 * as publication date) or label (required for values such as tissue that are backed by an ontology term ID).
 * @param key - Value to sort category value views by.
 * @returns Function that returns a number indicating sort precedence of cv0 vs cv1.
 */
function sortCategoryValueViews(key: "key" | "label") {
  return (
    cvv0: SelectCategoryValueView,
    cvv1: SelectCategoryValueView
  ): number => {
    return COLLATOR_CASE_INSENSITIVE.compare(cvv0[key], cvv1[key]);
  };
}

/**
 * Sort category views by display label, ascending.
 * @param c0 - First category views to compare.
 * @param c1 - Second category views to compare.
 * @returns Number indicating sort precedence of c0 vs c1.
 */
function sortCategoryViews(c0: CategoryView, c1: CategoryView): number {
  return COLLATOR_CASE_INSENSITIVE.compare(c0.label, c1.label);
}

/**
 * Summarize each category by applying the filters applicable to each category and counting the occurrences of category
 * values in each resulting result set.
 * @param originalRows - Original result set before filtering.
 * @param queries - Selected filters applicable to a category.
 * @returns Intermediate filter state with category, category values and counts fully built. Note, empty category values
 * and category value selected states are added after this initial structure is built.
 */
function summarizeCategories<T extends Categories>(
  originalRows: Row<T>[],
  queries: Query<T>[]
): FilterState {
  return queries.reduce((accum: FilterState, query: Query<T>) => {
    // Apply the filters on the original result set
    const rows = applyFilters(originalRows, query.filters);

    // Count the category value occurrences in each category that shares this filter. Range categories are not
    // summarized; they always use the full range of the original result set.
    query.categoryFilterIds.forEach((categoryFilterId: CATEGORY_FILTER_ID) => {
      if (!isMatchKindBetween(categoryFilterId)) {
        accum[categoryFilterId] = summarizeSelectCategory(
          categoryFilterId,
          rows
        );
      }
    });
    return accum;
  }, {} as FilterState);
}

/**
 * Count occurrences of category values across the result set for the given single or multiselect category.
 * @param categoryFilterId - ID of category to count category values of.
 * @param filteredRows - Array of rows containing category values to count.
 * @return Map of category values keyed by category value key.
 */
function summarizeSelectCategory<T extends Categories>(
  categoryFilterId: CATEGORY_FILTER_ID,
  filteredRows: Row<T>[]
): KeyedSelectCategoryValue {
  // Aggregate category value counts for each row.
  return filteredRows.reduce((accum: KeyedSelectCategoryValue, row: Row<T>) => {
    // Grab the values of the category for this dataset row.
    let categoryValues = row.values[categoryFilterId];

    if (!Array.isArray(categoryValues)) {
      categoryValues = [categoryValues];
    }

    // Init category value if it doesn't already exist. Default selected state to false (selected and selectedParital
    // state is updated from the filter state at a later point).
    categoryValues.forEach((categoryValueKey: CategoryValueId) => {
      let categoryValue = accum.get(categoryValueKey);
      if (!categoryValue) {
        categoryValue = {
          count: 0,
          key: categoryValueKey,
          selected: false,
          selectedPartial: false,
        };
        accum.set(categoryValueKey, categoryValue);
      }
      // Increment category value count.
      categoryValue.count++;
    });
    return accum;
  }, new Map<CategoryValueId, SelectCategoryValue>());
}

/**
 * Update category value as selected if it is not currently selected, otherwise remove the selected category value if
 * it's already selected.
 * @param selectedValue - ID of the selected category value.
 * @param currentSelectedCategoryValueKeys - Keys of the current set of selected category values.
 * @param multiselect - True if category allows more than one selected value.
 * @returns Array of selected category values.
 */
function toggleCategoryValueSelected(
  selectedValue: CategoryValueId,
  currentSelectedCategoryValueKeys: CategoryValueId[],
  multiselect: boolean
): CategoryValueId[] {
  // Convert to set for ease of lookup and lookup efficiency.
  const selectedCategoryValueKeySet = new Set(currentSelectedCategoryValueKeys);

  // Toggle selected state of each selected value.
  if (selectedCategoryValueKeySet.has(selectedValue)) {
    selectedCategoryValueKeySet.delete(selectedValue);
  } else {
    // If category only allows single selected value, clear all other values
    if (!multiselect) {
      selectedCategoryValueKeySet.clear();
    }
    selectedCategoryValueKeySet.add(selectedValue);
  }

  return [...selectedCategoryValueKeySet.values()];
}

/**
 * Track select of the given ontology category and category value.
 * @param config - Configuration model of selected category.
 * @param categoryValueKey - Selected category value key (e.g. "HsapDv:0000003").
 * @param filters - Current set of selected category values (values) or ranges keyed by category (id).
 */
function trackCuratedOntologyCategoryValueSelected<T extends Categories>(
  config: CuratedOntologyCategoryFilterConfig,
  categoryValueKey: CategoryValueId,
  filters: Filters<T>
) {
  const { analyticsEvent, categoryFilterId, mask } = config;

  // No tracking if event isn't specified on category config.
  if (!analyticsEvent) {
    return;
  }

  // Only track the select (and not deselect) of category value.
  const categoryFilters = new Set(
    getCategoryFilter(categoryFilterId, filters)?.value as CategoryValueId[]
  );
  if (!categoryFilters.has(categoryValueKey)) {
    // Grab the analytics event for this category.

    // Find the node for the selected value.
    const ontologySpeciesKey = getOntologySpeciesKey(categoryValueKey);
    const ontologyRootNodes = mask[ontologySpeciesKey];
    if (!ontologyRootNodes) {
      return; // Error state - ontology does not exist.
    }
    const selectedOntologyNode = findOntologyNodeById(
      ontologyRootNodes,
      categoryValueKey
    );
    if (!selectedOntologyNode) {
      return; // Error state - ontology node with given ID does not exist.
    }

    // Build up payload for tracking event and send.
    track(analyticsEvent, {
      label: selectedOntologyNode.label,
      ontologyTermId: selectedOntologyNode.ontology_term_id,
    });
  }
}

/**
 * Track select of the given select category and category value.
 * @param config - Configuration model of selected category.
 * @param categoryValueKey - Selected category value key (e.g. "10 3' v2").
 * @param filters - Current set of selected category values (values) or ranges keyed by category (id).
 */
function trackSelectCategoryValueSelected<T extends Categories>(
  config: CategoryFilterConfig,
  categoryValueKey: CategoryValueId,
  filters: Filters<T>
) {
  const { analyticsEvent, categoryFilterId } = config;

  // No tracking if event isn't specified on category config.
  if (!analyticsEvent) {
    return;
  }

  // Only track the select (and not deselect) of category value.
  const categoryFilters = new Set(
    getCategoryFilter(categoryFilterId, filters)?.value as CategoryValueId
  );
  if (!categoryFilters.has(categoryValueKey)) {
    // Build up payload for tracking event and send.
    const payload = categoryValueKey;
    track(analyticsEvent, { payload });
  }
}
