import { useCallback, useEffect, useState } from "react";
import { Filters, Row } from "react-table";
import { Ontology } from "src/common/entities";
import {
  buildCuratedOntologyCategoryView,
  onFilterCuratedOntologyCategory,
} from "src/common/hooks/useCategoryFilter/common/curatedOntologyUtils";
import {
  buildMultiPanelCategoryView,
  buildMultiPanelCurrentFilters,
  buildMultiPanelSelectedUIState,
  buildMultiPanelUIState,
  isMultiPanelCategoryFilterConfig,
  onFilterMultiPanelCategory,
} from "src/common/hooks/useCategoryFilter/common/multiPanelOntologyUtils";
import {
  getCategoryFilter,
  isCategoryValueIdSet,
} from "src/common/hooks/useCategoryFilter/common/utils";
import {
  CATEGORY_FILTER_CONFIGS_BY_ID,
  COLLATOR_CASE_INSENSITIVE,
} from "src/components/common/Filter/common/constants";
import {
  Categories,
  CATEGORY_FILTER_ID,
  CategoryFilter,
  CategoryFilterConfig,
  CategoryFilterPanelConfig,
  CategorySet,
  CategoryValueId,
  CategoryView,
  CuratedOntologyCategoryFilterConfig,
  FilterKey,
  FilterState,
  KeyedSelectCategoryValue,
  MultiPanelCategoryFilterUIState,
  MultiPanelSelectedUIState,
  MultiPanelUIState,
  ON_FILTER_SOURCE,
  OnFilterFn,
  ONTOLOGY_VIEW_KEY,
  OntologyNode,
  OntologyTermSet,
  Range,
  RangeCategory,
  SelectCategoryValue,
} from "src/components/common/Filter/common/entities";
import { listOntologyTreeIds } from "src/components/common/Filter/common/utils";
import {
  addRangeCategories,
  buildRangeCategoryView,
  isSelectedValueRange,
  onFilterRangeCategory,
} from "./common/rangeUtils";
import {
  buildSelectCategoryView,
  onFilterSelectCategory,
} from "./common/selectUtils";

/**
 * Shape of return value from this useFilter hook.
 */
export interface FilterInstance {
  categoryViews: CategoryView[];
  multiPanelSelectedUIState: MultiPanelSelectedUIState;
  onFilter: OnFilterFn;
}

/**
 * Selected filters applicable to a category; used when deriving category value counts from current set of filters.
 * Identical queries can be shared by categories to reduce the number of result set filtering.
 */
interface Query<T extends Categories> {
  categoryFilterIds: CATEGORY_FILTER_ID[];
  filters: Filters<T>;
}

/**
 * react-table function to call when updating set of selected filters.
 */
/* eslint-disable-next-line @typescript-eslint/no-explicit-any -- function type as per react-table's setFilter. */
type SetFilterFn = (columnId: string, updater: any) => void;

/**
 * Faceted filter functionality over dataset metadata. "or" between values, "and" across categories.
 * @param originalRows - Original result set before filtering.
 * @param categoryFilterIds - Set of category IDs to include for this filter instance.
 * @param filters - Current set of selected category values (values) keyed by category (id).
 * @param setFilter - Function to update set of selected values for a category.
 * @param initialMultiPanelSelectedUIState - Selected state of category to set as an initial state.
 * @returns Object containing filter accessor (view model of filter state) and filter mutator (function to modify react-
 * table's internal filter state).
 */
export function useCategoryFilter<T extends Categories>(
  originalRows: Row<T>[],
  categoryFilterIds: Set<CATEGORY_FILTER_ID>,
  filters: Filters<T>,
  setFilter: SetFilterFn,
  initialMultiPanelSelectedUIState: MultiPanelSelectedUIState
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
  const [multiPanelUIState, setMultiPanelUIState] =
    useState<Map<CATEGORY_FILTER_ID, MultiPanelCategoryFilterUIState>>();

  // Serializable selected state of multi-panel category filters; facilitates save of filter UI state to local storage.
  const [multiPanelSelectedUIState, setMultiPanelSelectedUIState] =
    useState<MultiPanelSelectedUIState>(initialMultiPanelSelectedUIState);

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

    setOntologyTermLabelsById(
      keyOntologyTermLabelsById(originalRows, categoryFilterIds)
    );
  }, [originalRows, categoryFilterIds, ontologyTermLabelsById]);

  // Build up UI hierarchies for each multi-panel category filter, used to facilitate easy calculation of selected
  // and partially selected states.
  useEffect(() => {
    // Only set multi-panel state if there are rows to parse category values from, only initialize once.
    if (!originalRows.length || multiPanelUIState) {
      return;
    }

    const uiState = buildMultiPanelUIState(
      originalRows,
      categoryFilterIds,
      initialMultiPanelSelectedUIState
    );
    setMultiPanelUIState(uiState);
  }, [
    categoryFilterIds,
    initialMultiPanelSelectedUIState,
    multiPanelUIState,
    originalRows,
  ]);

  // Build next filter state on change of filter.
  useEffect(() => {
    // Must have category set and multi-panel UI state before next filter state can be calculated.
    if (!categorySet || !multiPanelUIState) {
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
      selectedCategoryValue: CategoryValueId | Range,
      selectedLabel: string | Range,
      source: ON_FILTER_SOURCE = ON_FILTER_SOURCE.FILTER
    ) => {
      if (!categorySet || !ontologyTermLabelsById || !multiPanelUIState) {
        return; // Error state - category set, ontology labels map and multi-panel UI state should be set at this point.
      }

      // Grab the configuration model for the selected category.
      const config = CATEGORY_FILTER_CONFIGS_BY_ID[categoryFilterId];

      // Handle range categories.
      if (isSelectedValueRange(selectedCategoryValue)) {
        const nextFilters = onFilterRangeCategory(
          config,
          selectedCategoryValue
        );
        setFilter(categoryFilterId, nextFilters);
        return;
      }

      if (!selectedCategoryValue) {
        return; // Error state - selected value must be defined for select and ontology categories.
      }

      // Handle ontology categories.
      if (isCuratedOntologyCategoryFilterConfig(config)) {
        const nextFilters = onFilterCuratedOntologyCategory(
          config,
          selectedCategoryValue,
          filters,
          categorySet
        );
        setFilter(categoryFilterId, nextFilters);
        return;
      }

      // Handle multi-panel categories.
      if (isMultiPanelCategoryFilterConfig(config)) {
        const [nextMultiPanelUIState, nextFilters] = onFilterMultiPanelCategory(
          config,
          selectedCategoryValue,
          selectedLabel as string,
          source,
          multiPanelUIState
        );
        setMultiPanelUIState(nextMultiPanelUIState);
        setMultiPanelSelectedUIState(
          buildMultiPanelSelectedUIState(nextMultiPanelUIState)
        );
        setFilter(categoryFilterId, nextFilters);
        return;
      }

      // Handle single or multiselect categories.
      const nextFilters = onFilterSelectCategory(
        config,
        selectedCategoryValue,
        filters
      );
      setFilter(categoryFilterId, nextFilters);
    },
    [categorySet, filters, ontologyTermLabelsById, setFilter, multiPanelUIState]
  );

  return {
    categoryViews: buildCategoryViews(
      filterState,
      multiPanelUIState,
      ontologyTermLabelsById
    ),
    multiPanelSelectedUIState,
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
      (categoryValueId: CategoryValueId) => {
        if (!categoryValuesByKey.has(categoryValueId)) {
          categoryValuesByKey.set(categoryValueId, {
            categoryValueId,
            count: 0,
            selected: false,
            selectedPartial: false,
          });
        }
      }
    );
  }
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
        categoryOntologyIds = listOntologyTreeIds(config.source);
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
 * Build view-specific models from filter state, to facilitate easy rendering.
 * @param filterState - Categories, category value and their counts with the current filter applied.
 * @param multiPanelUIState - Current set of category values that the user has selected in multi-panel category filters.
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
      if (config.viewKind === "SELECT") {
        return buildSelectCategoryView(
          config,
          rangeOrSelectValue as KeyedSelectCategoryValue,
          filterState,
          ontologyTermLabelsById
        );
      }

      // Build view model for curated ontology categories.
      if (config.viewKind === "CURATED_ONTOLOGY") {
        return buildCuratedOntologyCategoryView(
          config as CuratedOntologyCategoryFilterConfig,
          rangeOrSelectValue as KeyedSelectCategoryValue,
          filterState,
          ontologyTermLabelsById
        );
      }

      // Build view model for range categories.
      if (config.viewKind === "RANGE") {
        return buildRangeCategoryView(
          config,
          rangeOrSelectValue as RangeCategory
        );
      }

      // Build view model for multi-panel categories.
      return buildMultiPanelCategoryView(
        config,
        rangeOrSelectValue as KeyedSelectCategoryValue,
        multiPanelUIState,
        ontologyTermLabelsById
      );
    })
    .filter((categoryView): categoryView is CategoryView => !!categoryView)
    .sort(sortCategoryViews);
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
    // Adding type assertion here: ontology values are backed by ontology arrays.
    // @ts-expect-error -- resolve mismatch between Categories and FilterKey.
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
 * @param categoryFilterIds - Set of category filter IDs to include for this filter instance.
 * @returns Map of ontology term labels keyed by ID.
 */
function keyOntologyTermLabelsById<T extends Categories>(
  originalRows: Row<T>[],
  categoryFilterIds: Set<CATEGORY_FILTER_ID>
): Map<string, string> {
  const labelsById = new Map<string, string>();

  // Collect the set of term sets across all category filters.
  const termSets = [...categoryFilterIds].reduce(
    (accum: OntologyTermSet[], categoryFilterId: CATEGORY_FILTER_ID) => {
      const config = CATEGORY_FILTER_CONFIGS_BY_ID[categoryFilterId];

      // Add source of curated ontology categories.
      if (config.valueSourceKind === "CURATED") {
        accum.push(config.source);
        return accum;
      }

      // Add source for each panel of multi-panel categories.
      if (config.viewKind === "MULTI_PANEL") {
        config.panels.forEach((panel: CategoryFilterPanelConfig) => {
          if (panel.sourceKind === "CURATED") {
            accum.push(panel.source);
          }
        });
      }
      return accum;
    },
    [] as OntologyTermSet[]
  );

  // Key curated term sets.
  termSets.forEach((termSet: OntologyTermSet) => {
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
          buildMultiPanelCurrentFilters(multiPanelUIState)
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
    categoryValues.forEach((categoryValueId: CategoryValueId) => {
      let categoryValue = accum.get(categoryValueId);
      if (!categoryValue) {
        categoryValue = {
          categoryValueId,
          count: 0,
          selected: false,
          selectedPartial: false,
        };
        accum.set(categoryValueId, categoryValue);
      }
      // Increment category value count.
      categoryValue.count++;
    });
    return accum;
  }, new Map<CategoryValueId, SelectCategoryValue>());
}
