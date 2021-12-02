// Display-optimized structure of category and corresponding category values and counts.
import { useCallback, useMemo, useState } from "react";
import { Filters, FilterValue, Row } from "react-table";
import {
  CategoryKey,
  CATEGORY_KEY,
  DatasetRow,
} from "src/components/common/Filter/common/entities";

// Metadata values grouped by metadata key.
export interface CategoryView {
  key: CATEGORY_KEY;
  values: CategoryValueView[];
}

// Selected category values in a category.
interface CategoryFilter {
  id: string;
  value: FilterValue;
}

// Set of all category values in the full result set, keyed by their corresponding category.
type CategorySet = { [K in CATEGORY_KEY]: Set<CategoryValueKey> };

// Metadata value, selected state and a set containing the IDs of the result set rows it's associated with. That is,
// dataset ID when filtering datasets, collection ID when filtering collections.
export interface CategoryValue {
  key: CategoryValueKey;
  associatedWith: Set<string>;
  selected: boolean;
}

// Category values to be used as keys when building filter functionality. For example, "Homo sapiens" or
// "10X 3' v2 sequencing".
export type CategoryValueKey = string;

// View model of metadata value, selected state and count.
export interface CategoryValueView {
  key: CategoryValueKey;
  count: number;
  selected: boolean;
}

// Shape of return value from this useFilter hook.
export interface FilterInstance {
  categories: CategoryView[];
  onFilter: OnFilterFn;
}

// State backing filter functionality and calculations. Converted to view model for display.
type FilterState = {
  [K in CATEGORY_KEY]: Map<CategoryValueKey, CategoryValue>;
};

// Filterable metadata keys where the type of the corresponding value is Ontology. Currently, that is all metadata
// keys except is_primary_data.
export type OntologyCategoryKey = keyof Omit<
  Record<CATEGORY_KEY, string>,
  CATEGORY_KEY.IS_PRIMARY_DATA
>;

// Filterable metadata keys where the type of the corresponding value is not Ontology. Currently, that is
// is_primary_data only.
export type NonOntologyCategoryKey = keyof Pick<
  Record<CATEGORY_KEY, string>,
  CATEGORY_KEY.IS_PRIMARY_DATA
>;

// Selected filters applicable to a category; used when deriving category value counts from current set of filters.
// Identical queries can be shared by categories to reduce the number of result set filtering.
interface Query {
  categoryKeys: CategoryKey[];
  filters: Filters<DatasetRow>;
}

// Function invoked when selected state of a category value is toggled.
export type OnFilterFn = (
  categoryKey: CategoryKey,
  categoryValueKey: CategoryValueKey
) => void;

// TODO(cc) check re-renders (5?)

/**
 * Faceted filter functionality over dataset metadata. "or" between values, "and" across categories.
 * @param originalRows - Original result set before filtering.
 * @param filters - Current set of selected category values (values) keyed by category (id).
 * @param setFilter - Function to update set of selected values for a category.
 * @param groupById - Field to group rows by when summarizing category values. For example, "collection_id" when
 * filtering collections.
 * @returns Object containing filter accessor (view model of filter state) and filter mutator (function to modify react-
 * table's internal filter state).
 */
export function useCategoryFilter(
  originalRows: Row<DatasetRow>[],
  filters: Filters<DatasetRow>,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any -- function type as per react-table's setFilter.
  setFilter: (columnId: string, updater: any) => void,
  groupById: string
): FilterInstance {
  // Complete set of categories and category values for the result set.
  const [categorySet, setCategorySet] = useState<CategorySet>();

  // Core filter state facilitating build of complete set of categories, category values and counts for a filtered
  // result set.
  const [filterState, setFilterState] = useState<FilterState>();

  // Set up original, full set of categories and their values.
  useMemo(() => {
    // Only build category set once on load.
    if (categorySet) {
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
      categorySet,
      groupById
    );
    setFilterState(nextFilterState);
  }, [categorySet, filters, groupById, originalRows]);

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
 * Add back any category values that have been filtered out, and set their values to 0 and maintain their selected
 * state.
 * @param nextFilterState - Filter state currently being built due to change in filter.
 * @param filters - Current set of selected category values (values) keyed by category (id).
 * @param categorySet - Original, unfiltered sets of category values keyed by their category.
 */
function addEmptyCategoryValues(
  nextFilterState: FilterState,
  filters: Filters<DatasetRow>,
  categorySet: CategorySet
) {
  // Check filter state for each category for missing (empty) category values.
  for (const [categoryKey, categoryValuesByCategoryValue] of Object.entries(
    nextFilterState
  )) {
    // Grab the expected set of category values.
    const allCategoryValueKeys = categorySet[categoryKey as CategoryKey];

    // Grab the filters for this category; used to set the selected state of any missing category values.
    const categoryFilter = getCategoryFilter(
      categoryKey as CategoryKey,
      filters
    );

    // If expected category value is missing from this category's category values, add it back in with a count of 0.
    // Check filters to see if category value is currently selected.
    [...allCategoryValueKeys.values()].forEach(
      (categoryValueKey: CategoryValueKey) => {
        if (!categoryValuesByCategoryValue.has(categoryValueKey)) {
          const selected = Boolean(
            categoryFilter && categoryFilter.value.includes(categoryValueKey)
          );
          categoryValuesByCategoryValue.set(categoryValueKey, {
            associatedWith: new Set(),
            key: categoryValueKey,
            selected,
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
function buildCategorySet(originalRows: Row<DatasetRow>[]): CategorySet {
  // Build up category values for each category
  return Object.values(CATEGORY_KEY).reduce(
    (accum: CategorySet, categoryKey: CategoryKey) => {
      // Check category value for this category, in every row.
      originalRows.forEach((originalRow: Row<DatasetRow>) => {
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
          count: categoryValue.associatedWith.size,
          key: categoryValue.key,
          selected: categoryValue.selected,
        }))
        .sort(sortCategoryValueViews);
      // Return completed view model of this category.
      return {
        key: categoryKey as CategoryKey,
        values: categoryValueViews,
      };
    })
    .sort(sortCategoryViews);
}

/**
 * Build categories, category values and counts for the updated filter. For each category, build up category values
 * counts by counting occurrences of category values across rows. Maintain selected category values state from filters.
 * Retain category values with 0 counts from given category set.
 * @param originalRows - Original result set before filtering.
 * @param filters - Current set of selected category values (values) keyed by category (id).
 * @param categorySet - Original, unfiltered sets of category values keyed by their category.
 * @param groupById - Field to group rows by when summarizing category values. For example, "collection_id" when
 * filtering collections.
 * @returns New filter state generated from the current set of selected category values.
 */
function buildNextFilterState(
  originalRows: Row<DatasetRow>[],
  filters: Filters<DatasetRow>,
  categorySet: CategorySet,
  groupById: string
): FilterState {
  // Remove empty categories from filter. react-table maintains an empty error for categories that previously had
  // selected values.
  const sanitizedFilters = filters.filter(
    (filter: CategoryFilter) => filter.value.length > 0
  );

  // Build set of filters that are applicable to each category.
  const queries = buildQueries(sanitizedFilters);

  // Build up base filter state of categories, category values and counts.
  const nextFilterState = summarizeCategories(originalRows, queries, groupById);

  // Update selected flag for the selected category values.
  flagSelectedCategoryValues(nextFilterState, filters);

  // Always display category values even if their count is 0; add back any category values that have been filtered out.
  addEmptyCategoryValues(nextFilterState, filters, categorySet);

  return nextFilterState;
}

/**
 * Update selected state of category values to match the current set of selected filters.
 * @param nextFilterState - Filter state being built on select of filter.
 * @param filters - Current set of selected category values (values) keyed by category (id).
 */
function flagSelectedCategoryValues(
  nextFilterState: FilterState,
  filters: Filters<DatasetRow>
) {
  // Build convenience set of selected keys; interpolated category key and category value key.
  const filterKeys = filters.reduce((accum, filter) => {
    filter.value.forEach((value: string) => {
      accum.add(`${filter.id}${value}`);
    });
    return accum;
  }, new Set<string>());

  Object.keys(nextFilterState).forEach((key: string) => {
    const categoryFilterState = nextFilterState[key as CategoryKey];
    [...categoryFilterState.keys()].forEach(
      (categoryValueKey: CategoryValueKey) => {
        const categoryValue = categoryFilterState.get(categoryValueKey);
        if (categoryValue) {
          categoryValue.selected = filterKeys.has(`${key}${categoryValueKey}`);
        }
      }
    );
  });
}

/**
 * Find and return the selected values for the given category.
 * @param categoryKey - Key of category to find selected filters of.
 * @param filters - Current set of selected category values (values) keyed by category (id).
 * @returns Array of filters
 */
function getCategoryFilter(
  categoryKey: CategoryKey,
  filters: Filters<DatasetRow>
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
function includesSome(
  originalRows: Row<DatasetRow>[],
  filters: CategoryFilter[]
): Row<DatasetRow>[] {
  // Return all rows if there are no filters.
  if (filters.length === 0) {
    return originalRows;
  }
  return originalRows.filter((row: Row<DatasetRow>) => {
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
 * Determine the set of filters that are applicable to each category. That is, for a category, all selected filters
 * other than the selected filters for that category can be applied to the result set to determine the counts for
 * @param filters - Current set of selected category values (values) keyed by category (id).
 * @returns Array of query models representing of the selected filters applicable for each category.
 */
function buildQueries(filters: Filters<DatasetRow>): Query[] {
  return Object.values(CATEGORY_KEY).reduce(
    (accum: Query[], categoryKey: CategoryKey) => {
      // Determine the filters that are applicable to this category.
      const filtersExcludingSelf = filters.filter((filter: CategoryFilter) => {
        return filter.id !== categoryKey;
      });

      // Check if we have an existing  query with an identical filter. If so, add category to that query. Otherwise
      // create new query for this filter.
      const matchingQuery = accum.find((query: Query) =>
        isFilterEqual(query.filters, filtersExcludingSelf)
      );
      if (matchingQuery) {
        (matchingQuery as Query).categoryKeys.push(categoryKey);
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
 * Determine if given filters are identical.
 * @param filters0 - First filter to compare.
 * @param filters1 - Second filter to compare.
 */
function isFilterEqual(
  filters0: Filters<DatasetRow>,
  filters1: Filters<DatasetRow>
): boolean {
  return (
    filters0.length === filters1.length &&
    filters0.every((val) => filters1.includes(val))
  );
}

/**
 * Build updated set of selected filters for the given category and the selected category value.
 * @param categoryKey - Category key (i.e. "disease") of selected category value.
 * @param categoryValueKey - Category value key (e.g. "normal") to toggle selected state of.
 * @param filters - Current set of selected category values.
 * @returns Array of selected category values for the given category.
 */
function buildNextCategoryFilters(
  categoryKey: CategoryKey,
  categoryValueKey: CategoryValueKey,
  filters: Filters<DatasetRow>
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
 * Sort category views by key, ascending.
 * @param c0 - First category to compare.
 * @param c1 - Second category to compare.
 * @returns Number indicating sort precedence of c0 vs c1.
 */
function sortCategoryViews(c0: CategoryView, c1: CategoryView): number {
  if (c0.key < c1.key) {
    return -1;
  }
  if (c0.key > c1.key) {
    return 1;
  }
  return 0;
}

/**
 * Sort category value views by key, ascending.
 * @param cvv0 - First filtered rows to compare.
 * @param cvv1 - Second filtered rows to compare.
 * @returns Number indicating sort precedence of cv0 vs cv1.
 */
function sortCategoryValueViews(
  cvv0: CategoryValueView,
  cvv1: CategoryValueView
): number {
  if (cvv0.key < cvv1.key) {
    return -1;
  }
  if (cvv0.key > cvv1.key) {
    return 1;
  }
  return 0;
}

/**
 * Count occurrences of category values across the result set for the given category. Count is modelled as the set of
 * "group by" values on each category value. For datasets, the grouping is per dataset (so dataset ID). For collections,
 * the grouping is per collection (so collection_id).
 * @param categoryKey - Category to count category values.
 * @param filteredRows - Array of rows containing category values to count.
 * @param groupById - Field to group rows by when summarizing category values. For example, "collection_id" when
 * filtering collections.
 * @return Map of category values keyed by category value key.
 */
function summarizeCategory(
  categoryKey: CategoryKey,
  filteredRows: Row<DatasetRow>[],
  groupById: string
): Map<CategoryValueKey, CategoryValue> {
  // Aggregate category value counts for each row.
  return filteredRows.reduce(
    (accum: Map<CategoryValueKey, CategoryValue>, row: Row<DatasetRow>) => {
      // Grab the value of the category for this dataset row, convert to array if it isn't already (for example,
      // is_primary_data).
      const rawCategoryValue = row.values[categoryKey];
      const categoryValueKeys = Array.isArray(rawCategoryValue)
        ? rawCategoryValue
        : [rawCategoryValue];

      // Init category value it doesn't already exist. Default selected state to false (selected state is updated
      // from the filter state at a later point).
      categoryValueKeys.forEach((categoryValueKey: CategoryValueKey) => {
        let categoryValue = accum.get(categoryValueKey);
        if (!categoryValue) {
          categoryValue = {
            associatedWith: new Set<string>(),
            key: categoryValueKey,
            selected: false,
          };
          accum.set(categoryValueKey, categoryValue);
        }
        // Add group by value for this category value key.
        categoryValue.associatedWith.add(row.values[groupById]);
      });
      return accum;
    },
    new Map<CategoryValueKey, CategoryValue>()
  );
}

/**
 * Summarize each category by applying the filters applicable to each category and counting the occurrences of category
 * values in each resulting result set.
 * @param originalRows - Original result set before filtering.
 * @param queries - Selected filters applicable to a category.
 * @param groupById - Field to group rows by when summarizing category values. For example, "collection_id" when
 * filtering collections.
 * @returns Intermediate filter state with category, category values and counts fully built. Note, empty category values
 * and category value selected states are added after this initial structure is built.
 */
function summarizeCategories(
  originalRows: Row<DatasetRow>[],
  queries: Query[],
  groupById: string
): FilterState {
  return queries.reduce((accum: FilterState, query: Query) => {
    // Apply the filters on the original result set
    const rows = includesSome(originalRows, query.filters);

    // Count the category value occurrences in each category that shares this filter.
    query.categoryKeys.forEach((categoryKey: CategoryKey) => {
      accum[categoryKey] = summarizeCategory(categoryKey, rows, groupById);
    });

    return accum;
  }, {} as FilterState);
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
