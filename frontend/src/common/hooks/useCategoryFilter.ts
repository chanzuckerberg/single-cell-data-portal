// Display-optimized structure of category and corresponding category values and counts.
import { useCallback, useEffect, useState } from "react";
import { Filters, FilterValue, Row } from "react-table";
import { Ontology } from "src/common/entities";
import {
  CATEGORY_FILTER_CONFIGS_BY_ID,
  CATEGORY_FILTER_UI_CONFIGS,
  COLLATOR_CASE_INSENSITIVE,
  DEVELOPMENT_STAGE_ONTOLOGY_TERM_SET,
  TISSUE_ORGAN_ONTOLOGY_TERM_SET,
  TISSUE_SYSTEM_ONTOLOGY_TERM_SET,
} from "src/components/common/Filter/common/constants";
import {
  Categories,
  CategoryFilterConfig,
  CategoryFilterPanelConfig,
  CategoryValueKey,
  CategoryView,
  CategoryViews,
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
  OntologyTermSet,
  ONTOLOGY_VIEW_KEY,
  ONTOLOGY_VIEW_LABEL,
  OrFilterPrefix,
  ORGANISM,
  PUBLICATION_DATE_LABELS,
  Range,
  RangeCategoryView,
  SelectCategoryValueView,
  SelectCategoryView,
} from "src/components/common/Filter/common/entities";
import {
  findOntologyNodeById,
  findOntologyParentNode,
  getOntologySpeciesKey,
  listOntologyTreeIds,
} from "src/components/common/Filter/common/utils";
import { track } from "../analytics";

/**
 * Entry in react-table's filters arrays, models selected category values in a category.
 */
interface CategoryFilter {
  id: string;
  value: FilterValue;
}

/**
 * Filterable metadata object key. For example, "assay" or "cell_type". Used for object key lookups.
 * TODO(cc) move to entities or delete even? rename usages (ie variable names) to match
 */
export type CategoryFilterId = keyof Record<CATEGORY_FILTER_ID, string>;

/*
 * Set of all category values in the full result set, keyed by their corresponding category.
 */
type CategorySet = { [K in CATEGORY_FILTER_ID]: CategorySetValue };

/**
 * Possible category set values, either a set of category key values (for single or multiselect categories, or ontology
 * categories) or a range.
 */
type CategorySetValue = Set<CategoryValueKey> | Range;

/**
 * Internal filter model of a single or multiselect category value, or an ontology category value: category value keyed
 * by category value key (for easy look-up when summarizing category).
 */
type KeyedSelectCategoryValue = Map<CategoryValueKey, SelectCategoryValue>;

/**
 * Internal filter model of a single or multiselect category, or an ontology category.
 */
export interface SelectCategoryValue {
  key: CategoryValueKey;
  count: number;
  selected: boolean;
}

/**
 * Shape of return value from this useFilter hook.
 */
export interface FilterInstance {
  categoryViews: CategoryViews[];
  onFilter: OnFilterFn;
}

/**
 * State backing filter functionality and calculations. Converted to view model for display.
 */
type FilterState = {
  [K in CATEGORY_FILTER_ID]: RangeCategory | KeyedSelectCategoryValue;
};

/**
 * Selected filters applicable to a category; used when deriving category value counts from current set of filters.
 * Identical queries can be shared by categories to reduce the number of result set filtering.
 */
interface Query<T extends Categories> {
  categoryFilterIds: CategoryFilterId[];
  filters: Filters<T>;
}

/**
 * Internal filter model of a range category.
 */
interface RangeCategory {
  key: CategoryValueKey;
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
      categorySet
    );
    setFilterState(nextFilterState);
  }, [categoryFilterIds, categorySet, filters, originalRows]);

  // Update set of filters on select of category value. Track selected category value.
  const onFilter = useCallback<OnFilterFn>(
    (
      categoryFilterId: CategoryFilterId,
      selectedValue: CategoryValueKey | Range
    ) => {
      if (!categorySet) {
        return; // Error state - category set should be set at this point.
      }

      // Grab the configuration model for the selected category.
      const config = CATEGORY_FILTER_CONFIGS_BY_ID[categoryFilterId];

      // Handle range categories.
      if (!isCategoryValueKey(selectedValue)) {
        onFilterRangeCategory(config, selectedValue, setFilter);
        return;
      }

      // Handle ontology categories.
      if (isCuratedOntologyCategoryFilterConfig(config)) {
        onFilterOntologyCategory(
          config,
          selectedValue,
          setFilter,
          filters,
          categorySet
        );
        return;
      }

      // Handle single or multiselect categories.
      const nextCategoryFilters = onFilterSelectCategory(
        config,
        selectedValue,
        filters
      );
      setFilter(categoryFilterId, nextCategoryFilters);
    },
    [categorySet, filters, setFilter]
  );

  return {
    categoryViews: buildCategoryViews(filterState, ontologyTermLabelsById),
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
      categorySet[categoryFilterId as CategoryFilterId];
    if (
      !allCategoryValueKeys || // Error state - all category values for this category can't be found.
      !isCategorySetCategoryKeyValue(allCategoryValueKeys) // Error state - should be category key value.
    ) {
      return;
    }

    // If expected category value is missing from this category's category values, add it back in with a count of 0.
    [...allCategoryValueKeys.values()].forEach(
      (categoryValueKey: CategoryValueKey) => {
        if (!categoryValuesByKey.has(categoryValueKey)) {
          categoryValuesByKey.set(categoryValueKey, {
            count: 0,
            key: categoryValueKey,
            selected: false,
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
  selectedCategoryValues: Set<CategoryValueKey>,
  categoryKeyValues: Set<CategoryValueKey>
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
      isCategorySetCategoryKeyValue(categorySetRange) // Error state - should be a range.
    ) {
      return;
    }

    // Add range to next filter state.
    const { min: originalMin, max: originalMax } = categorySetRange;
    nextFilterState[categoryKey] = {
      key: categoryKey,
      max: originalMax ?? 0,
      min: originalMin ?? 0,
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
      if (isMatchKindBetween(filter.id as CategoryFilterId)) {
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
    (accum: CategorySet, categoryFilterId: CategoryFilterId) => {
      // Calculate the initial state of range categories.
      if (isMatchKindBetween(categoryFilterId)) {
        const counts = originalRows
          // Use filterKey to pull value from original row
          .map((originalRow) => originalRow.values[categoryFilterId])
          .filter((count) => !!count || count === 0); // Remove bad data, just in case!

        accum[categoryFilterId] = {
          max: counts.length ? Math.max(...counts) : 0,
          min: counts.length ? Math.min(...counts) : 0,
        };
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
        let categoryValueSet = accum[categoryFilterId] as Set<CategoryValueKey>;
        if (!categoryValueSet) {
          categoryValueSet = new Set<CategoryValueKey>();
          accum[categoryFilterId] = categoryValueSet;
        }
        // Add the category values for this row to the set.
        let values: CategoryValueKey | CategoryValueKey[] =
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
          .filter((value: CategoryValueKey) => {
            // If category is an ontology, confirm value is included in ontology tree for display.
            return !isCategoryOntology || categoryOntologyIds.has(value);
          })
          .forEach((value: CategoryValueKey) => categoryValueSet.add(value));
      });
      return accum;
    },
    {} as CategorySet
  );
}

/**
 * Build the display value for the given category and category value.
 * @param config - Config model of category to build category value views for.
 * @param categoryValueKey - Category value to display (e.g. "normal").
 * @param ontologyTermLabelsById - Set of ontology term labels keyed by term ID, used to determine labels for ontology
 * TODO(cc) docs, drilling
 * @returns String to display as a label for the given category and category value.
 */
function buildCategoryValueLabel(
  config: CategoryFilterConfig,
  categoryValueKey: CategoryValueKey,
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

  // TODO(cc) revisit this - needs drilling of config and ontology term labels
  if (config.labelKind === "LOOKUP_LABEL_BY_TERM_ID") {
    // De-tag values. TODO(cc) add config for this. revisit logic.
    let processedCategoryValueKey = categoryValueKey;
    if (config.categoryFilterId === "TISSUE_CALCULATED") {
      const tokens = categoryValueKey.split(/:(.*)/s);
      processedCategoryValueKey = tokens[1];
    }

    return (
      ontologyTermLabelsById.get(processedCategoryValueKey) ?? categoryValueKey
    ); // TODO(cc) error handling here?
  }

  // Return all other category values as is.
  return categoryValueKey;
}

/**
 * Build view-specific models from UI config and filter state, to facilitate easy rendering.
 * @param filterState - Categories, category value and their counts with the current filter applied.
 * @param ontologyTermLabelsById - Set of ontology term labels keyed by term ID, used to determine labels for ontology
 * terms.
 * @returns Array of category views objects.
 */
function buildCategoryViews(
  filterState?: FilterState,
  ontologyTermLabelsById?: Map<string, string>
): CategoryViews[] {
  if (!filterState || !ontologyTermLabelsById) {
    return [];
  }

  // Build up view models for each UI config.
  // TODO(cc) pass in original set of category IDs and build UI from that (rather than having to do truthy check for categoryValueByValue below)
  return CATEGORY_FILTER_UI_CONFIGS.map((categoryFilterUIConfig) => {
    const { categoryFilterConfigIds, label } = categoryFilterUIConfig;
    const categoryViews = categoryFilterConfigIds
      .map((categoryConfigKey) => {
        const config = CATEGORY_FILTER_CONFIGS_BY_ID[categoryConfigKey];
        const { categoryFilterId } = config;

        // Return if there's no filter state for the UI config. This will be true for categories such as cell count and
        // gene count when viewing collections.
        const rangeOrSelectValue = filterState[categoryFilterId];
        if (!rangeOrSelectValue) {
          return;
        }

        // Build view model for single or multiselect categories
        // TODO(cc) remove need for multiple checks, can we get just kind to assert the category value type?
        if (
          config.viewKind === "SELECT" &&
          isSelectCategoryValue(rangeOrSelectValue)
        ) {
          return buildSelectCategoryView(
            categoryFilterId, // TODO(cc) remove id from params, here and others
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
          return buildOntologyCategoryView(
            categoryFilterId,
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
          return buildRangeCategoryView(
            categoryFilterId,
            config,
            rangeOrSelectValue
          );
        }

        // Build view model for multi-panel categories.
        // TODO(cc)
        if (
          config.viewKind === "MULTI_PANEL" &&
          isSelectCategoryValue(rangeOrSelectValue)
        ) {
          return buildMultiPanelCategoryView(
            categoryFilterId,
            config,
            rangeOrSelectValue,
            ontologyTermLabelsById
          );
        }

        // Error - unknown view kind/category value.
        console.log(
          `Error attempting to build view model for category "${config.categoryFilterId}".`
        );
      })
      .filter((categoryView): categoryView is CategoryView => !!categoryView);

    return {
      categoryViews,
      label,
    };
  }).sort(sortCategoryViews);
}

// /**
//  * Restrict related categories where applicable. For example, tissue system restricts selectable values in
//  * tissue organ and tissue.
//  * @param nextFilterState - Filter state currently being built due to change in filter.
//  * @param filters - Current set of selected category values (values) keyed by category (id).
//  */
// function applyCrossCategoryRestrictions<T extends Categories>(
//   nextFilterState: FilterState,
//   filters: Filters<T>
// ) {
//   return nextFilterState;
//   // // Apply restrictions to categories that are restricted by another category.
//   // // TODO(cc) generalize and move to function
//   // Object.keys(nextFilterState).forEach((categoryFilterId: string) => {
//   //   const config =
//   //     CATEGORY_FILTER_CONFIGS_BY_ID[categoryFilterId as CategoryFilterId];
//   //
//   //   // Ignore category filters that are not restricted by a parent.
//   //   if (config.valueRestrictionKind !== "CHILDREN_OF_SELECTED_PARENT_TERMS") {
//   //     return;
//   //   }
//   //
//   //   // If there are no selected values in parents, then there are no restrictions to be applied.
//   //   const isAnyParentValueSelected = config.parentCategoryFilterIds.some(
//   //     (parentCategoryFilterId) =>
//   //       !!getCategoryFilter(parentCategoryFilterId, filters)
//   //   );
//   //   if (!isAnyParentValueSelected) {
//   //     return;
//   //   }
//   //
//   //   // Iterate through parent category filters and determine the set of allowed values.
//   //   const includeValues = new Set();
//   //   for (const parentCategoryFilterId of config.parentCategoryFilterIds) {
//   //     // If there are no selected values for this category, check the next parent.
//   //     const categoryFilter = getCategoryFilter(parentCategoryFilterId, filters);
//   //     if (!categoryFilter) {
//   //       continue;
//   //     }
//   //
//   //     // Determine the set of allowed values; add the descendants for every selected value in parent category filter,
//   //     // unless the selected value has a descendant that is already in the allowed values. If we add both the parent
//   //     // selected value (in addition to a descendant) then the allowed set of values will over-show.
//   //     categoryFilter.value.forEach((selectedValue: string) => {
//   //       const descendants = DESCENDANTS[selectedValue] ?? [];
//   //       const isAnyDescendantSelected = descendants.some((descendant) =>
//   //         includeValues.has(descendant)
//   //       );
//   //       if (isAnyDescendantSelected) {
//   //         return;
//   //       }
//   //
//   //       // No descendant is already added to the set of allowed values; add all descendants as allowed values.
//   //       descendants.forEach((descendant) => {
//   //         includeValues.add(descendant);
//   //       });
//   //     });
//   //   }
//   //
//   //   // Remove any values that are not in the allowed set.
//   //   const categoryFilterState =
//   //     nextFilterState[categoryFilterId as CategoryFilterId];
//   //
//   //   // TODO(cc) fix
//   //   if (!isSelectCategoryValue(categoryFilterState)) {
//   //     return;
//   //   }
//   //
//   //   [...categoryFilterState.keys()].forEach(
//   //     (categoryValueKey: CategoryValueKey) => {
//   //       if (
//   //         !includeValues.has(categoryValueKey) &&
//   //         !categoryFilterState.get(categoryValueKey)?.selected // TODO(cc) revisit optional chaining
//   //       ) {
//   //         categoryFilterState.delete(categoryValueKey);
//   //       }
//   //     }
//   //   );
//   // });
// }

/**
 * Build categories, category values and counts for the updated filter. For each category, build up category values
 * counts by counting occurrences of category values across rows. Maintain selected category values state from filters.
 * Retain category values with 0 counts from given category set.
 * @param originalRows - Original result set before filtering.
 * @param categoryFilterIds - Set of category IDs to include for this filter instance.
 * @param filters - Current set of selected category values (values) keyed by category (id).
 * @param categorySet - Original, unfiltered sets of category values keyed by their category.
 * @returns New filter state generated from the current set of selected category values.
 */
function buildNextFilterState<T extends Categories>(
  originalRows: Row<T>[],
  categoryFilterIds: Set<CATEGORY_FILTER_ID>,
  filters: Filters<T>,
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
  setSelectedStates(nextFilterState, filters);

  // Restrict related categories where applicable. For example, tissue system restricts selectable values in
  // tissue organ and tissue.
  // TODO(cc) remove?
  // applyCrossCategoryRestrictions(nextFilterState, filters);

  return nextFilterState;
}

/**
 * Build updated set of selected filters for the given ontology tree category and the selected category value.
 * @param categoryFilterId - Category ID (i.e. "development stage") of selected category value.
 * @param categoryValueKey - Selected category value key (e.g. "HsapDv:0000003") to update selected state of.
 * @param filters - Current set of selected category values.
 * @param categoryKeyValues - Original, full set of values for this category.
 * @param mask - View model of ontology for this category.
 * @returns Array of selected category values for the given category.
 */
export function buildNextOntologyCategoryFilters<T extends Categories>(
  categoryFilterId: CategoryFilterId,
  categoryValueKey: CategoryValueKey,
  filters: Filters<T>,
  categoryKeyValues: Set<CategoryValueKey>,
  mask: OntologyTermSet
): CategoryValueKey[] {
  // Grab the current selected values for the category.
  const categoryFilters = new Set(
    getCategoryFilter(categoryFilterId, filters)?.value as CategoryValueKey[]
  );

  // TODO(cc) resolve with #2569 - will need to change OntologyView key or add "translate category value key" functionality to getOntologySpeciesKey (to handle "CL:xxx" for tissue)
  // TODO(cc) do we still need this?
  const speciesKey = /* categoryFilterId === CATEGORY_FILTER_ID.TISSUE
      ? ONTOLOGY_VIEW_KEY.UBERON
      :*/ categoryValueKey;

  // Find the selected and parent node, if any, for the selected value.
  const ontologySpeciesKey = getOntologySpeciesKey(speciesKey);
  const ontologyRootNodes = mask[ontologySpeciesKey];
  if (!ontologyRootNodes) {
    return [...categoryFilters.values()]; // Error state - ontology does not exist.
  }
  const selectedOntologyNode = findOntologyNodeById(
    ontologyRootNodes,
    categoryValueKey
  );
  if (!selectedOntologyNode) {
    return [...categoryFilters.values()]; // Error state - ontology node with given ID does not exist.
  }
  const parentNode = findOntologyParentNode(
    ontologyRootNodes,
    selectedOntologyNode
  );

  // Toggle selected state of selected category value.
  if (categoryFilters.has(categoryValueKey)) {
    // Selected value is already in the set of selected values, remove it.
    categoryFilters.delete(categoryValueKey);

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
    categoryFilters.add(categoryValueKey);

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
 * Build updated set of selected filters for the given single or multiselect category and the selected category value.
 * Ontology categories are handled separately.
 * @param categoryFilterId - ID (i.e. "disease") of selected category value.
 * @param categoryValueKey - Category value key (e.g. "normal") to toggle selected state of.
 * @param filters - Current set of selected category values.
 * @returns Array of selected category values for the given category.
 * TODO(cc) revert?
 */
function buildNextSelectCategoryFilters<T extends Categories>(
  categoryFilterId: CategoryFilterId,
  categoryValueKey: CategoryValueKey,
  filters: Filters<T>
): CategoryValueKey[] {
  // Grab the current selected values for the category.
  const categoryFilters = getCategoryFilter(categoryFilterId, filters);

  // Currently, no filters already selected for this category; add category value as first.
  if (!categoryFilters) {
    return [categoryValueKey];
  }

  // Create new array of selected category value keys, with the selected state of the given category value toggled.
  const multiselect =
    CATEGORY_FILTER_CONFIGS_BY_ID[categoryFilterId].multiselect;
  return toggleCategoryValueSelected(
    categoryValueKey,
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
    (accum: Query<T>[], categoryFilterId: CategoryFilterId) => {
      // Determine the filters that are applicable to this category.
      const config = CATEGORY_FILTER_CONFIGS_BY_ID[categoryFilterId];
      const applicableFilters = filters.filter((filter: CategoryFilter) => {
        // Handle standard faceted search where filters are not filtered by their own selected values.
        if (config.queryKind === "EXCLUDES_SELF") {
          return filter.id !== categoryFilterId;
        }

        // TODO(cc) remove (and remove kind)?
        // // Handle category filters that are neither restricted by themselves nor restricted by their children.
        // const { childrenCategoryFilterIds } = config;
        // return (
        //   filter.id !== categoryFilterId &&
        //   !childrenCategoryFilterIds.includes(filter.id as CATEGORY_FILTER_ID)
        // );
      });

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
 * Build view model of ontology category.
 * @param categoryFilterId - ID of category to find selected filters of.
 * @param config - Config model of a curated ontology category.
 * @param categoryValueByValue - Internal filter model of ontology category.
 * @param filterState - Categories, category value and their counts with the current filter applied. Required when
 * checking enabled state of view that is dependent on the state of another category.
 * @returns Ontology view model.
 * TODO(cc) revisit naming - this is specific to curated ontologies
 * TODO(cc) docs, drilling
 */
function buildOntologyCategoryView(
  categoryFilterId: CategoryFilterId,
  config: CuratedOntologyCategoryFilterConfig,
  categoryValueByValue: KeyedSelectCategoryValue,
  filterState: FilterState,
  ontologyTermLabelsById: Map<string, string>
): OntologyCategoryView {
  const { isLabelVisible, isSearchable, isZerosVisible, label, mask } = config;

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
 * Build view model of node of ontology tree to be displayed as a value in an ontology menu.
 * @param ontologyNode - Ontology node to build view model for.
 * @param categoryValueByValue - Internal filter model of ontology category.
 * @returns Ontology view model.
 * TODO(cc) docs, drilling
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
    };
  }

  // Build up base view model.
  const view = {
    count: categoryValue.count,
    key: categoryValueKey,
    label: ontologyTermLabelsById.get(categoryValueKey) ?? categoryValueKey,
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

/**
 * Build view model of single or multiselect category.
 * @param categoryFilterId - ID of category to find selected filters of.
 * @param config - Config model of ontology category.
 * @param categoryValueByValue - Internal filter model of single or multiselect category.
 * checking enabled state of view that is dependent on the state of another category.
 * @param ontologyTermLabelsById - Set of ontology term labels keyed by term ID, used to determine labels for ontology
 * @returns Select category view model.
 */
function buildMultiPanelCategoryView(
  categoryFilterId: CategoryFilterId,
  config: OntologyMultiPanelFilterConfig,
  categoryValueByValue: KeyedSelectCategoryValue,
  ontologyTermLabelsById: Map<string, string> // TODO(cc) revisit drilling
): OntologyMultiPanelCategoryView {
  // Build select view models for all values to be displayed. TODO(cc) only build for I and then don't need dedupe?
  const allSelectCategoryValueViews = buildSelectCategoryValueViews(
    config,
    [...categoryValueByValue.values()],
    ontologyTermLabelsById
  ).sort(sortCategoryValueViews);

  // // Remove duplicated view models created from removal of inferred and explicit prefixes on filter values. Always
  // // take the larger count when deduping view models as the larger count will include both inferred and explicit values.
  // const selectCategoryValueViews = dedupePrefixedSelectCategoryValueViews(
  //   allSelectCategoryValueViews
  // );

  // Build up view model for each panel specified in config.
  const { panels } = config;
  const panelViews = panels.map((panel: CategoryFilterPanelConfig) => {
    return {
      label: panel.label,
      views:
        panel.sourceKind === "ONLY_CURATED"
          ? maskOnlyCurated(panel.mask, allSelectCategoryValueViews)
          : maskExcludingCurated(
              panel.excludeMasks,
              allSelectCategoryValueViews
            ),
    };
  });

  // Build view model of multi-panel category.
  const multiPanelView: OntologyMultiPanelCategoryView = {
    key: categoryFilterId,
    label: config.label,
    panels: panelViews,
  };

  // TODO(cc) pinned values
  // TODO(cc) check disabled state here

  return multiPanelView;
}

/**
 * Remove curated and left-over ancestors.
 * TODO(cc) docs, location, name
 */
function maskExcludingCurated(
  masks: OntologyTermSet[],
  selectCategoryValueViews: SelectCategoryValueView[]
): SelectCategoryValueView[] {
  const excludeValues = masks.reduce((accum, mask) => {
    accum.push(...[...listOntologyTreeIds(mask)]);
    return accum;
  }, [] as string[]);
  return selectCategoryValueViews.filter((view) => {
    const [prefix, processedKey] = view.key.split(/:(.*)/s); // TODO(cc) utils, revisit prefix
    if (excludeValues.includes(processedKey)) {
      return false;
    }
    return prefix !== OrFilterPrefix.INFERRED;
  });
}

/**
 * TODO(cc) docs, location, name
 */
function maskOnlyCurated(
  mask: OntologyTermSet,
  selectCategoryValueViews: SelectCategoryValueView[]
): SelectCategoryValueView[] {
  const allowedValues = [...listOntologyTreeIds(mask)];
  return selectCategoryValueViews.filter((view) => {
    const [, processedKey] = view.key.split(/:(.*)/s); // TODO(cc) utils
    return allowedValues.includes(processedKey);
  });
}

/**
 * Remove duplicated view models created from removal of inferred and explicit prefixes on filter values. Always
 * take the larger count when deduping view models as the larger count will include both inferred and explicit values.
 * TODO(cc) only required for tissue and cell type, should we configure this? or just do it for all? docs plus location.
 */
function dedupePrefixedSelectCategoryValueViews(
  selectCategoryValueViews: SelectCategoryValueView[]
): SelectCategoryValueView[] {
  const dedupedSelectCateogryValueViews = selectCategoryValueViews.reduce(
    (accum, selectCategoryValueView) => {
      const { key, count } = selectCategoryValueView;
      const [, processedKey] = key.split(/:(.*)/s); // TODO(cc) utils
      if (accum.has(processedKey)) {
        // eslint-disable-next-line @typescript-eslint/no-non-null-assertion -- using accum.has(key) above to ensure accum has value.
        const addedSelectCategoryValueView = accum.get(processedKey)!;
        const { count: currentCount } = addedSelectCategoryValueView;
        const max = Math.max(currentCount, count);
        if (max !== currentCount) {
          addedSelectCategoryValueView.count = max;
        }
      } else {
        accum.set(processedKey, selectCategoryValueView);
      }
      return accum;
    },
    new Map<string, SelectCategoryValueView>()
  );

  return [...dedupedSelectCateogryValueViews.values()];
}

/**
 * Build view model of range category.
 * @param categoryFilterId - ID of category to find selected filters of.
 * @param config - Config model of ontology category.
 * @param rangeCategory - Internal filter model of range category.
 * @returns Range view model.
 */
function buildRangeCategoryView(
  categoryFilterId: CategoryFilterId,
  config: CategoryFilterConfig,
  rangeCategory: RangeCategory
): RangeCategoryView {
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
 * Build select category value view models of the given selected values.
 * @param config - Config model of category to build category value views for.
 * @param categoryValues - Values to build view models from.
 * @param ontologyTermLabelsById - Set of ontology term labels keyed by term ID, used to determine labels for ontology
 */
function buildSelectCategoryValueViews(
  config: CategoryFilterConfig,
  categoryValues: SelectCategoryValue[],
  ontologyTermLabelsById: Map<string, string> // TODO(cc) revisit drilling, maybe have a format fn for label instead?
): SelectCategoryValueView[] {
  return categoryValues.map(({ count, key, selected }: SelectCategoryValue) => {
    return {
      count,
      key,
      label: buildCategoryValueLabel(config, key, ontologyTermLabelsById),
      selected: selected,
    };
  });
}

/**
 * Build view model of single or multiselect category.
 * @param categoryFilterId - ID of category to find selected filters of.
 * @param config - Config model of ontology category.
 * @param categoryValueByValue - Internal filter model of single or multiselect category.
 * @param filterState - Categories, category value and their counts with the current filter applied. Required when
 * @param ontologyTermLabelsById - Set of ontology term labels keyed by term ID, used to determine labels for ontology
 * checking enabled state of view that is dependent on the state of another category.
 * @returns Select category view model.
 */
function buildSelectCategoryView(
  categoryFilterId: CategoryFilterId,
  config: CategoryFilterConfig,
  categoryValueByValue: KeyedSelectCategoryValue,
  filterState: FilterState,
  ontologyTermLabelsById: Map<string, string> // TODO(cc) revisit drilling
): SelectCategoryView {
  // Grab the config for this category.
  const { pinnedCategoryValues, tooltip } =
    CATEGORY_FILTER_CONFIGS_BY_ID[categoryFilterId];

  // Remove any excluded values from the filter.
  // TODO(cc) revisit location of this. move earlier? update comments here but essentially the prevents an organ from appearing in the "catch all" tissue.
  // let selectCategoryValues = [...categoryValueByValue.values()];
  // if (config.valueRestrictionKind === "CHILDREN_OF_SELECTED_PARENT_TERMS") {
  //   const excludeTerms = config.parentCategoryFilterIds.reduce(
  //     (accum: string[], parentCategoryFilterId) => {
  //       const parentConfig =
  //         CATEGORY_FILTER_CONFIGS_BY_ID[parentCategoryFilterId];
  //       if (parentConfig.valueSourceKind !== "CURATED") {
  //         return accum; // Error state - parent must be mask? TODO(cc) revisit does it have to be a mask? or can it be another field?
  //       }
  //       accum.push(...listOntologyTreeIds(parentConfig.mask));
  //       return accum;
  //     },
  //     []
  //   );
  //
  //   if (excludeTerms && excludeTerms.length) {
  //     selectCategoryValues = selectCategoryValues.filter(
  //       (selectCategoryValue) =>
  //         !selectCategoryValue.key ||
  //         !excludeTerms.includes(selectCategoryValue.key)
  //     );
  //   }
  // }

  const allCategoryValueViews = buildSelectCategoryValueViews(
    config,
    [...categoryValueByValue.values()],
    ontologyTermLabelsById
  ).sort(sortCategoryValueViews);

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
  categoryFilterId: CategoryFilterId,
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
  selectedCategoryValues: Set<CategoryValueKey>,
  categoryKeyValues: Set<CategoryValueKey>
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
  selectedCategoryValues: Set<CategoryValueKey>
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
 * TODO(cc) revisit - this was never type narrowing? can we improve this? also remove config truthy check?
 */
function isMatchKindBetween(categoryFilterId: CategoryFilterId): boolean {
  const config = CATEGORY_FILTER_CONFIGS_BY_ID[categoryFilterId];
  if (!config) {
    console.log(`Config not found for category "${categoryFilterId}".`);
    return false; // Error state - return false.
  }
  return config.matchKind === "BETWEEN";
}

/**
 * Determine if the given ethnicity is considered unspecified (that is, na or unknown).
 * @param categoryValueKey - Ethnicity value to check if it's specified.
 * @returns True if ethnicity is either na or unknown.
 */

function isEthnicitySpecified(categoryValueKey: CategoryValueKey) {
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
  selectedCategoryValues: Set<CategoryValueKey>,
  categoryKeyValues: Set<CategoryValueKey>
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
 * Determine if the given selected value is a selected category value key (and not a range).
 * @param selectedValue - Selected filter value, either a category value key (e.g. "normal") or a range (for example,
 * {min: n, max: n}).
 * @returns True if given selected value is a selected category value.
 * TODO(cc) revisit
 */
function isCategoryValueKey(
  selectedValue: CategoryValueKey | Range
): selectedValue is CategoryValueKey {
  return !(
    (selectedValue as Range).min >= 0 && (selectedValue as Range).max >= 0
  );
}

/**
 * Determine if the given category value is a select category value (and not a range category value).
 * @param categorySetValue - Range or category set value.
 * @returns True if category set value is a set of category value keys.
 * TODO(cc) revisit
 */
function isCategorySetCategoryKeyValue(
  categorySetValue: CategorySetValue
): categorySetValue is Set<CategoryValueKey> {
  return categorySetValue instanceof Set;
}

/**
 * Determine if the given category config is an ontology (and not a "regular" category config).
 * @param config - Config model of category, either an ontology category config or a base category config.
 * TODO(cc) revisit - naming, use etc
 */
function isCuratedOntologyCategoryFilterConfig(
  config: CategoryFilterConfig
): config is CuratedOntologyCategoryFilterConfig {
  return config.valueSourceKind == "CURATED";
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
 * determine if development stage species should be visible.
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
 * TODO(cc) iterate through all nodes to add term ID/term label to map
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
 * TODO(cc) key all values in curated ontology
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
 * TODO(cc) key row ontology term values
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
  pinnedCategoryValues?: CategoryValueKey[]
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
  selectedCategoryValues: Set<CategoryValueKey>
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
 * Handle select of ontology value: build and set next set of filters for this category. Track selected ontology value.
 * @param config - Configuration model of selected category.
 * @param selectedValue - Selected category value key (e.g. [1, 100]).
 * @param setFilter - Function to update set of selected values for a category.
 * @param filters - Current set of selected category values (values) or ranges keyed by category (id).
 * @param categorySet - Original, unfiltered sets of category values keyed by their category.
 */
function onFilterOntologyCategory<T extends Categories>(
  config: CuratedOntologyCategoryFilterConfig, // TODO(cc) revisit this, should it be any ontology type?
  selectedValue: CategoryValueKey,
  setFilter: SetFilterFn,
  filters: Filters<T>,
  categorySet: CategorySet
) {
  const { categoryFilterId, mask } = config;

  // Track selected category and value.
  trackOntologyCategoryValueSelected(config, selectedValue, filters);

  // Build and set next set of filters for this category.
  const nextCategoryFilters = buildNextOntologyCategoryFilters(
    categoryFilterId,
    selectedValue,
    filters,
    categorySet[categoryFilterId] as Set<CategoryValueKey>,
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
  const { max, min } = selectedValue;
  console.log(max, min);
  // Track select of new range mim/max, ignoring any clear of selected range. Only track if event is specified on
  // configuration model
  if (analyticsEvent /* && selectedValue.length > 0*/) {
    track(analyticsEvent, {
      max,
      min,
    });
  }

  // Update filters for this range category. Convert range min/max object to tuple for react-table.
  setFilter(categoryFilterId, [min, max]);
}

/**
 * Handle select of select value: build and set next set of filters for this category. Track selected select value.
 * @param config - Configuration model of selected category.
 * @param selectedValue - Selected category value key (e.g. [1, 100]).
 * @param setFilter - Function to update set of selected values for a category.
 * @param filters - Current set of selected category values (values) or ranges keyed by category (id).
 * TODO(cc) docs, change others to match
 */
function onFilterSelectCategory<T extends Categories>(
  config: CategoryFilterConfig,
  selectedValue: CategoryValueKey,
  // setFilter: SetFilterFn,
  filters: Filters<T>
): string[] {
  const { categoryFilterId } = config;

  // Track selected category and value.
  trackSelectCategoryValueSelected(config, selectedValue, filters);

  // Build and set next set of filters for this category.
  return buildNextSelectCategoryFilters(
    categoryFilterId,
    selectedValue,
    filters
  );
}

/**
 * Update selected state of categories to match the current set of selected filters.
 * @param nextFilterState - Filter state being built on select of filter.
 * @param filters - Current set of selected category values (values) or ranges keyed by category (id).
 */
function setSelectedStates<T extends Categories>(
  nextFilterState: FilterState,
  filters: Filters<T>
) {
  Object.keys(nextFilterState).forEach((categoryKey: string) => {
    // Grab the filter state for this category.
    const categoryFilterState =
      nextFilterState[categoryKey as CategoryFilterId];

    // Grab the filters for this category.
    const categoryFilter = getCategoryFilter(
      categoryKey as CategoryFilterId,
      filters
    );
    if (!categoryFilter || !categoryFilter.value) {
      // There are no selected values for this category
      return;
    }

    // Handle single and multiselect categories, or ontology categories.
    if (isSelectCategoryValue(categoryFilterState)) {
      // Create set for easy lookup of category values.
      const selectedCategoryValues = new Set(categoryFilter.value);

      // Check each category value in this category to see if it's selected.
      for (const [
        categoryValueKey,
        categoryValue,
      ] of categoryFilterState.entries()) {
        categoryValue.selected = selectedCategoryValues.has(categoryValueKey);
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
 * Sort category value views by key, ascending.
 * @param cvv0 - First category value view to compare.
 * @param cvv1 - Second category value view to compare.
 * @returns Number indicating sort precedence of cv0 vs cv1.
 */
function sortCategoryValueViews(
  cvv0: SelectCategoryValueView,
  cvv1: SelectCategoryValueView
): number {
  return COLLATOR_CASE_INSENSITIVE.compare(cvv0.label, cvv1.label);
}

/**
 * Sort category views by display label, ascending.
 * @param c0 - First category views to compare.
 * @param c1 - Second category views to compare.
 * @returns Number indicating sort precedence of c0 vs c1.
 */
function sortCategoryViews(c0: CategoryViews, c1: CategoryViews): number {
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
    query.categoryFilterIds.forEach((categoryFilterId: CategoryFilterId) => {
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
  categoryFilterId: CategoryFilterId,
  filteredRows: Row<T>[]
): KeyedSelectCategoryValue {
  // Aggregate category value counts for each row.
  return filteredRows.reduce((accum: KeyedSelectCategoryValue, row: Row<T>) => {
    // Grab the values of the category for this dataset row.
    let categoryValues = row.values[categoryFilterId];

    if (!Array.isArray(categoryValues)) {
      categoryValues = [categoryValues];
    }

    // Init category value if it doesn't already exist. Default selected state to false (selected state is updated
    // from the filter state at a later point).
    categoryValues.forEach((categoryValueKey: CategoryValueKey) => {
      let categoryValue = accum.get(categoryValueKey);
      if (!categoryValue) {
        categoryValue = {
          count: 0,
          key: categoryValueKey,
          selected: false,
        };
        accum.set(categoryValueKey, categoryValue);
      }
      // Increment category value count.
      categoryValue.count++;
    });
    return accum;
  }, new Map<CategoryValueKey, SelectCategoryValue>());
}

/**
 * Update category value as selected if it is not currently selected, otherwise remove the selected category value if
 * it's already selected.
 * @param selectedCategoryValueKey - Key of the selected category.
 * @param selectedCategoryValueKeys - Keys of the current set of selected category values.
 * @param multiselect - True if category allows more than one selected value.
 * @returns Array of selected category values.
 */
function toggleCategoryValueSelected(
  selectedCategoryValueKey: CategoryValueKey,
  selectedCategoryValueKeys: CategoryValueKey[],
  multiselect: boolean
): CategoryValueKey[] {
  // Convert to set for ease of lookup and lookup efficiency.
  const selectedCategoryValueKeySet = new Set(selectedCategoryValueKeys);
  if (selectedCategoryValueKeySet.has(selectedCategoryValueKey)) {
    selectedCategoryValueKeySet.delete(selectedCategoryValueKey);
  } else {
    // If category only allows single selected value, clear all other values
    if (!multiselect) {
      selectedCategoryValueKeySet.clear();
    }
    selectedCategoryValueKeySet.add(selectedCategoryValueKey);
  }
  return [...selectedCategoryValueKeySet.values()];
}

/**
 * Track select of the given ontology category and category value.
 * @param config - Configuration model of selected category.
 * @param categoryValueKey - Selected category value key (e.g. "HsapDv:0000003").
 * @param filters - Current set of selected category values (values) or ranges keyed by category (id).
 */
function trackOntologyCategoryValueSelected<T extends Categories>(
  config: CuratedOntologyCategoryFilterConfig, // TODO(cc) - revisit type here
  categoryValueKey: CategoryValueKey,
  filters: Filters<T>
) {
  const { analyticsEvent, categoryFilterId, mask } = config;

  // No tracking if event isn't specified on category config.
  if (!analyticsEvent) {
    return;
  }

  // Only track the select (and not deselect) of category value.
  const categoryFilters = new Set(
    getCategoryFilter(categoryFilterId, filters)?.value as CategoryValueKey[]
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
  categoryValueKey: CategoryValueKey,
  filters: Filters<T>
) {
  const { analyticsEvent, categoryFilterId } = config;

  // No tracking if event isn't specified on category config.
  if (!analyticsEvent) {
    return;
  }

  // Only track the select (and not deselect) of category value.
  const categoryFilters = new Set(
    getCategoryFilter(categoryFilterId, filters)?.value as CategoryValueKey[]
  );
  if (!categoryFilters.has(categoryValueKey)) {
    // Build up payload for tracking event and send.
    const payload = categoryValueKey;
    track(analyticsEvent, { payload });
  }
}
