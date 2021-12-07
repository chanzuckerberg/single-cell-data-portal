// Display-optimized structure of category and corresponding category values and counts.
import { useCallback, useMemo, useState } from "react";
import { Filters, FilterValue, Row } from "react-table";
import { COLLATOR_CASE_INSENSITIVE } from "src/components/common/Filter/common/constants";
import {
  Categories,
  CategoryValueKey,
  CategoryValueView,
  CategoryView,
  CATEGORY_KEY,
  CATEGORY_LABEL,
  IS_PRIMARY_DATA_LABEL,
  OnFilterFn,
} from "src/components/common/Filter/common/entities";

/* Entry in react-table's filters arrays, models selected category values in a category. */
interface CategoryFilter {
  id: string;
  value: FilterValue;
}

/* Filterable metadata object key. For example, "assay", "cell_type" or "is_primary_data". Used for object key lookups */
export type CategoryKey = keyof Record<CATEGORY_KEY, string>;

/* Set of all category values in the full result set, keyed by their corresponding category. */
type CategorySet = { [K in CATEGORY_KEY]: Set<CategoryValueKey> };

/* Metadata value, count and selected state. */
export interface CategoryValue {
  key: CategoryValueKey;
  count: number;
  selected: boolean;
}

/* Shape of return value from this useFilter hook. */
export interface FilterInstance {
  categories: CategoryView[];
  onFilter: OnFilterFn;
}

/* State backing filter functionality and calculations. Converted to view model for display. */
type FilterState = {
  [K in CATEGORY_KEY]: Map<CategoryValueKey, CategoryValue>;
};

/* Selected filters applicable to a category; used when deriving category value counts from current set of filters.
   Identical queries can be shared by categories to reduce the number of result set filtering. */
interface Query<T extends Categories> {
  categoryKeys: CategoryKey[];
  filters: Filters<T>;
}

/* react-table function to call when updating set of selected filters. */
/* eslint-disable-next-line @typescript-eslint/no-explicit-any -- function type as per react-table's setFilter. */
type SetFilterFn = (columnId: string, updater: any) => void;

/**
 * Faceted filter functionality over dataset metadata. "or" between values, "and" across categories.
 * @param originalRows - Original result set before filtering.
 * @param filters - Current set of selected category values (values) keyed by category (id).
 * @param setFilter - Function to update set of selected values for a category.
 * @returns Object containing filter accessor (view model of filter state) and filter mutator (function to modify react-
 * table's internal filter state).
 */
export function useCategoryFilter<T extends Categories>(
  originalRows: Row<T>[],
  filters: Filters<T>,
  setFilter: SetFilterFn
): FilterInstance {
  // Complete set of categories and category values for the result set.
  const [categorySet, setCategorySet] = useState<CategorySet>();

  // Core filter state facilitating build of complete set of categories, category values and counts for a filtered
  // result set.
  const [filterState, setFilterState] = useState<FilterState>();

  // Set up original, full set of categories and their values.
  useMemo(() => {
    // Only build category set if there are rows to parse category values from. Only build category set once on load.
    if (!originalRows.length || categorySet) {
      return;
    }
    setCategorySet(buildCategorySet(originalRows));
  }, [originalRows, categorySet]);

  // Build next filter state on change of filter.
  useMemo(() => {
    // Must have category set before next filter state can be calculated.
    if (!categorySet) {
      return;
    }
    const nextFilterState = buildNextFilterState(
      originalRows,
      filters,
      categorySet
    );
    setFilterState(nextFilterState);
  }, [categorySet, filters, originalRows]);

  // Update set of filters on select of category value.
  const onFilter = useCallback<OnFilterFn>(
    (categoryKey: CategoryKey, categoryValueKey: CategoryValueKey) => {
      const nextCategoryFilters = buildNextCategoryFilters(
        categoryKey,
        categoryValueKey,
        filters
      );
      setFilter(categoryKey, nextCategoryFilters);
    },
    [filters, setFilter]
  );

  return {
    categories: buildCategoryViews(filterState),
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
  for (const [categoryKey, categoryValuesByKey] of Object.entries(
    nextFilterState
  )) {
    // Grab the expected set of category values.
    const allCategoryValueKeys = categorySet[categoryKey as CategoryKey];
    if (!allCategoryValueKeys) {
      return; // Error state - all category values for this category can't be found.
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
 * Set up model of original, complete set of categories and their values.
 * @param originalRows - Original result set before filtering.
 * @returns Sets of category values keyed by their category.
 */
function buildCategorySet<T extends Categories>(
  originalRows: Row<T>[]
): CategorySet {
  // Build up category values for each category
  return Object.values(CATEGORY_KEY).reduce(
    (accum: CategorySet, categoryKey: CategoryKey) => {
      // Check category value for this category, in every row.
      originalRows.forEach((originalRow: Row<T>) => {
        // Grab the category values already added for this category, create new set if it hasn't already been created.
        let categoryValueSet = accum[categoryKey];
        if (!categoryValueSet) {
          categoryValueSet = new Set<CategoryValueKey>();
          accum[categoryKey] = categoryValueSet;
        }
        // Add the category values for this row to the set.
        let values: CategoryValueKey | CategoryValueKey[] =
          originalRow.values[categoryKey];
        if (!values) {
          console.log(`No values found for category "${categoryKey}".`);
          return accum;
        }
        if (!Array.isArray(values)) {
          values = [values];
        }
        values.forEach((value: CategoryValueKey) =>
          categoryValueSet.add(value)
        );
      });
      return accum;
    },
    {} as CategorySet
  );
}

/**
 * Build the display value for the given category and category value.
 * @param categoryKey - Category of category value (e.g. "disease").
 * @param categoryValueKey - Category value to display (e.g. "normal").
 * @returns String to display as a label for the given category and category value.
 */
function buildCategoryValueLabel(
  categoryKey: CategoryKey,
  categoryValueKey: CategoryValueKey
): string {
  // Only is_primary_data category values are to be transformed.
  if (categoryKey !== CATEGORY_KEY.IS_PRIMARY_DATA) {
    return categoryValueKey;
  }

  return IS_PRIMARY_DATA_LABEL[
    categoryValueKey as keyof typeof IS_PRIMARY_DATA_LABEL
  ];
}

/**
 * Build view-specific models from filter state, to facilitate easy rendering.
 * @param filterState - Categories, category value and their counts with the current filter applied.
 * @returns Array of category view object.s
 */
function buildCategoryViews(filterState?: FilterState): CategoryView[] {
  if (!filterState) {
    return [];
  }
  return Object.keys(filterState)
    .map((categoryKey: string) => {
      // Build category value view models for this category and sort.
      const categoryValueByValue = filterState[categoryKey as CategoryKey];
      const categoryValueViews = [...categoryValueByValue.values()]
        .map((categoryValue: CategoryValue) => ({
          count: categoryValue.count,
          key: categoryValue.key,
          label: buildCategoryValueLabel(
            categoryKey as CategoryKey,
            categoryValue.key
          ),
          selected: categoryValue.selected,
        }))
        .sort(sortCategoryValueViews);
      // Return completed view model of this category.
      return {
        key: categoryKey as CategoryKey,
        label: CATEGORY_LABEL[categoryKey as CategoryKey],
        values: categoryValueViews,
      };
    })
    .sort(sortCategoryViews);
}

/**
 * Build updated set of selected filters for the given category and the selected category value.
 * @param categoryKey - Category key (i.e. "disease") of selected category value.
 * @param categoryValueKey - Category value key (e.g. "normal") to toggle selected state of.
 * @param filters - Current set of selected category values.
 * @returns Array of selected category values for the given category.
 */
function buildNextCategoryFilters<T extends Categories>(
  categoryKey: CategoryKey,
  categoryValueKey: CategoryValueKey,
  filters: Filters<T>
): CategoryValueKey[] {
  // Grab the current selected values for the category.
  const categoryFilters = getCategoryFilter(categoryKey, filters);

  // Current no filters already selected for this category; add category value as first.
  if (!categoryFilters) {
    return [categoryValueKey];
  }

  // Create new array of selected category value keys, with the selected state of the given category value toggled.
  return toggleCategoryValueSelected(categoryValueKey, categoryFilters.value);
}

/**
 * Build categories, category values and counts for the updated filter. For each category, build up category values
 * counts by counting occurrences of category values across rows. Maintain selected category values state from filters.
 * Retain category values with 0 counts from given category set.
 * @param originalRows - Original result set before filtering.
 * @param filters - Current set of selected category values (values) keyed by category (id).
 * @param categorySet - Original, unfiltered sets of category values keyed by their category.
 * @returns New filter state generated from the current set of selected category values.
 */
function buildNextFilterState<T extends Categories>(
  originalRows: Row<T>[],
  filters: Filters<T>,
  categorySet: CategorySet
): FilterState {
  // Build set of filters that are applicable to each category.
  const queries = buildQueries(filters);

  // Build up base filter state of categories, category values and counts.
  const nextFilterState = summarizeCategories(originalRows, queries);

  // Always display category values even if their count is 0; add back any category values that have been filtered out.
  addEmptyCategoryValues(nextFilterState, categorySet);

  // Update selected flag for the selected category values.
  flagSelectedCategoryValues(nextFilterState, filters);

  return nextFilterState;
}

/**
 * Determine the set of filters that are applicable to each category. That is, for a category, all selected filters
 * other than the selected filters for that category can be applied to the result set to determine the counts for.
 * @param filters - Current set of selected category values (values) keyed by category (id).
 * @returns Array of query models representing of the selected filters applicable for each category.
 */
function buildQueries<T extends Categories>(filters: Filters<T>): Query<T>[] {
  return Object.values(CATEGORY_KEY).reduce(
    (accum: Query<T>[], categoryKey: CategoryKey) => {
      // Determine the filters that are applicable to this category.
      const filtersExcludingSelf = filters.filter((filter: CategoryFilter) => {
        return filter.id !== categoryKey;
      });

      // Check if we have an existing  query with an identical filter. If so, add category to that query. Otherwise
      // create new query for this filter.
      const matchingQuery = accum.find((query: Query<T>) =>
        isFilterEqual(query.filters, filtersExcludingSelf)
      );
      if (matchingQuery) {
        (matchingQuery as Query<T>).categoryKeys.push(categoryKey);
      } else {
        accum.push({
          categoryKeys: [categoryKey],
          filters: filtersExcludingSelf,
        });
      }
      return accum;
    },
    []
  );
}

/**
 * Update selected state of category values to match the current set of selected filters.
 * @param nextFilterState - Filter state being built on select of filter.
 * @param filters - Current set of selected category values (values) keyed by category (id).
 */
function flagSelectedCategoryValues<T extends Categories>(
  nextFilterState: FilterState,
  filters: Filters<T>
) {
  Object.keys(nextFilterState).forEach((categoryKey: string) => {
    // Grab the filter state for this category.
    const categoryFilterState = nextFilterState[categoryKey as CategoryKey];

    // Grab the filters for this category.
    const categoryFilter = getCategoryFilter(
      categoryKey as CategoryKey,
      filters
    );
    if (!categoryFilter || !categoryFilter.value) {
      // There are no selected values for this category
      return;
    }

    // Create set for easy lookup of category values.
    const selectedCategoryValues = new Set(categoryFilter.value);

    // Check each category value in this category to see if it's selected.
    for (const [
      categoryValueKey,
      categoryValue,
    ] of categoryFilterState.entries()) {
      categoryValue.selected = selectedCategoryValues.has(categoryValueKey);
    }
  });
}

/**
 * Find and return the selected values for the given category.
 * @param categoryKey - Key of category to find selected filters of.
 * @param filters - Current set of selected category values (values) keyed by category (id).
 * @returns Array of filters
 */
function getCategoryFilter<T extends Categories>(
  categoryKey: CategoryKey,
  filters: Filters<T>
): CategoryFilter | undefined {
  return filters.find((filter) => filter.id === categoryKey);
}

/**
 * Determine the rows that have values matching the given filters. Row must have at least one selected value across each
 * category. Mimics react-query's includeSome functionality.
 * @param originalRows - Original result set before filtering.
 * @param filters - Selected filters to apply to rows.
 * @returns Filtered array of rows.
 */
function includesSome<T extends Categories>(
  originalRows: Row<T>[],
  filters: CategoryFilter[]
): Row<T>[] {
  // Return all rows if there are no filters.
  if (filters.length === 0) {
    return originalRows;
  }
  return originalRows.filter((row: Row<T>) => {
    // "and" across categories.
    return filters.every((filter: CategoryFilter) => {
      const rowValue = row.values[filter.id];
      // "or" across category values (that is, inside a category).
      return (
        rowValue &&
        rowValue.length &&
        filter.value.some((val: string) => rowValue.includes(val)) // Handles string or array values
      );
    });
  });
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
 * Sort category value views by key, ascending.
 * @param cvv0 - First category value view to compare.
 * @param cvv1 - Second category value view to compare.
 * @returns Number indicating sort precedence of cv0 vs cv1.
 */
function sortCategoryValueViews(
  cvv0: CategoryValueView,
  cvv1: CategoryValueView
): number {
  return COLLATOR_CASE_INSENSITIVE.compare(cvv0.key, cvv1.key);
}

/**
 * Sort category views by key, ascending.
 * @param c0 - First category view to compare.
 * @param c1 - Second category view to compare.
 * @returns Number indicating sort precedence of c0 vs c1.
 */
function sortCategoryViews(c0: CategoryView, c1: CategoryView): number {
  return COLLATOR_CASE_INSENSITIVE.compare(c0.key, c1.key);
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
    const rows = includesSome(originalRows, query.filters);

    // Count the category value occurrences in each category that shares this filter.
    query.categoryKeys.forEach((categoryKey: CategoryKey) => {
      accum[categoryKey] = summarizeCategory(categoryKey, rows);
    });

    return accum;
  }, {} as FilterState);
}

/**
 * Count occurrences of category values across the result set for the given category.
 * @param categoryKey - Category to count category values.
 * @param filteredRows - Array of rows containing category values to count.
 * @return Map of category values keyed by category value key.
 */
function summarizeCategory<T extends Categories>(
  categoryKey: CategoryKey,
  filteredRows: Row<T>[]
): Map<CategoryValueKey, CategoryValue> {
  // Aggregate category value counts for each row.
  return filteredRows.reduce(
    (accum: Map<CategoryValueKey, CategoryValue>, row: Row<T>) => {
      // Grab the values of the category for this dataset row.
      const categoryValues = row.values[categoryKey];

      // Init category value it doesn't already exist. Default selected state to false (selected state is updated
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
    },
    new Map<CategoryValueKey, CategoryValue>()
  );
}

/**
 * Update category value as selected if it is not currently selected, otherwise remove the selected category value if
 * it's already selected.
 * @param selectedCategoryValueKey - Key of the selected category.
 * @param selectedCategoryValueKeys - Keys of the current set of selected category values.
 * @returns Array of selected category values.
 */
function toggleCategoryValueSelected(
  selectedCategoryValueKey: CategoryValueKey,
  selectedCategoryValueKeys: CategoryValueKey[]
): CategoryValueKey[] {
  // Convert to set for ease of lookup and lookup efficiency.
  const selectedCategoryValueKeySet = new Set(selectedCategoryValueKeys);
  if (selectedCategoryValueKeySet.has(selectedCategoryValueKey)) {
    selectedCategoryValueKeySet.delete(selectedCategoryValueKey);
  } else {
    selectedCategoryValueKeySet.add(selectedCategoryValueKey);
  }
  return [...selectedCategoryValueKeySet.values()];
}
