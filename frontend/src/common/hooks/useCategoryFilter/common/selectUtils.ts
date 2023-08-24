import { Filters } from "react-table";
import { track } from "src/common/analytics";
import {
  buildCategoryValueLabel,
  getCategoryFilter,
  sortCategoryValueViews,
} from "src/common/hooks/useCategoryFilter/common/utils";
import {
  CATEGORY_FILTER_CONFIGS_BY_ID,
  TOOLTIP_CATEGORY_DISABLED,
} from "src/components/common/Filter/common/constants";
import {
  Categories,
  CategoryFilterConfig,
  CategoryValueId,
  CATEGORY_FILTER_ID,
  FilterState,
  KeyedSelectCategoryValue,
  ORGANISM,
  SelectCategoryValue,
  SelectCategoryValueView,
  SelectCategoryView,
} from "src/components/common/Filter/common/entities";

/**
 * Utils specific to select category filters.
 */

/**
 * Build view model of single or multiselect category.
 * @param config - Config model of ontology category.
 * @param categoryValueByValue - Internal filter model of single or multiselect category.
 * @param filterState - Categories, category value and their counts with the current filter applied. Required when
 * @param ontologyTermLabelsById - Set of ontology term labels keyed by term ID, used to determine labels for ontology
 * checking enabled state of view that is dependent on the state of another category.
 * @returns Select category view model.
 */
export function buildSelectCategoryView(
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
  ).sort(sortCategoryValueViews("categoryValueId"));

  // Split values into pinned and unpinned.
  const [pinnedValues, unpinnedValues] = partitionSelectCategoryValueViews(
    allCategoryValueViews,
    pinnedCategoryValues
  );

  // Build view model of select category.
  const selectView: SelectCategoryView = {
    categoryFilterId: categoryFilterId,
    label: config.label,
    pinnedValues,
    unpinnedValues,
    values: allCategoryValueViews,
  };

  // Handle special cases where select category may be disabled.
  if (
    categoryFilterId === CATEGORY_FILTER_ID.SELF_REPORTED_ETHNICITY &&
    !isSelfReportedEthnicityViewEnabled(filterState)
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
 * Build updated set of selected filters for the given single or multiselect category and the selected category values.
 * @param config - Configuration model of selected category.
 * @param selectedValue - Category value ID to toggle the selected state of.
 * @param currentFilters - Current set of selected category values.
 * @returns Array of selected category values for the given category.
 */
export function buildNextSelectCategoryFilters<T extends Categories>(
  config: CategoryFilterConfig,
  selectedValue: CategoryValueId,
  currentFilters: Filters<T>
): CategoryValueId[] {
  const { categoryFilterId, multiselect } = config;

  // Grab the current selected values for the category.
  const categoryFilters = getCategoryFilter(categoryFilterId, currentFilters);

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
 * Build select category value view models of the given select category values.
 * @param config - Config model of category to build category value views for.
 * @param selectCategoryValues - Values to build view models from.
 * @param ontologyTermLabelsById - Set of ontology term labels keyed by term ID, used to determine labels for ontology
 * @param visibleUINodeIds - an optional set of visibleNode ids used when creating SelectCategoryValueViews for the multi panel ontology filters.
 * @returns an array of SelectCategoryValueViews for a category filter.
 */
export function buildSelectCategoryValueViews(
  config: CategoryFilterConfig,
  selectCategoryValues: SelectCategoryValue[],
  ontologyTermLabelsById: Map<string, string>,
  visibleUINodeIds?: Set<string>
): SelectCategoryValueView[] {
  return selectCategoryValues.map(
    ({
      categoryValueId,
      count,
      selected,
      selectedPartial,
    }: SelectCategoryValue) => {
      return {
        categoryValueId,
        count,
        label: buildCategoryValueLabel(
          config,
          categoryValueId,
          ontologyTermLabelsById
        ),
        selected,
        selectedPartial,
        visible: visibleUINodeIds
          ? visibleUINodeIds.has(categoryValueId)
          : true,
      };
    }
  );
}

/**
 * The self-reported ethnicity filter is only view enabled if:
 * 1. Homo sapiens is available as an option in the organism filter.
 * 2. The organism filter has selected values that includes Homo sapiens.
 * @param filterState - Categories, category value and their counts with the current filter applied. Required to
 * determine if self-reported ethnicity category should be enabled.
 * @returns True if self-reported ethnicity is either na or unknown.
 */
function isSelfReportedEthnicityViewEnabled(filterState: FilterState) {
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
    .map((selectCategoryValue) => selectCategoryValue.categoryValueId);
  return (
    selectedOrganisms.length === 0 ||
    selectedOrganisms.includes(ORGANISM.HOMO_SAPIENS)
  );
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
 * Handle select of select value: build and set next set of filters for this category. Track selected select value.
 * @param config - Configuration model of selected category.
 * @param selectedValue - Category value ID to use as selected value.
 * @param filters - Current set of selected category values (values) or ranges keyed by category (id).
 * @returns The array of selected values with latest selected applied.
 */
export function onFilterSelectCategory<T extends Categories>(
  config: CategoryFilterConfig,
  selectedValue: CategoryValueId,
  filters: Filters<T>
) {
  // Track selected category and value.
  trackSelectCategoryValueSelected(config, selectedValue, filters);

  // Build and set next set of filters for this category.
  return buildNextSelectCategoryFilters(config, selectedValue, filters);
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
    SelectCategoryValueView[],
  ] = [[], []];
  return categoryValues.reduce((accum, categoryValue) => {
    const [pinned, nonPinned] = accum;
    if (pinnedCategoryValues.includes(categoryValue.categoryValueId)) {
      pinned.push(categoryValue);
    } else {
      nonPinned.push(categoryValue);
    }
    return accum;
  }, partitionedValues);
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
  const { analyticsEvent, analyticsPayloadKey, categoryFilterId } = config;

  // No tracking if event isn't specified on category config.
  if (!analyticsEvent) {
    return;
  }

  // Only track the select (and not deselect) of category value.
  const categoryFilters = new Set(
    getCategoryFilter(categoryFilterId, filters)?.value as CategoryValueId
  );
  if (!categoryFilters.has(categoryValueKey)) {
    const payload = buildAnalyticsPayload(
      categoryValueKey,
      analyticsPayloadKey
    );
    track(analyticsEvent, payload);
  }
}

/**
 * Build up payload for tracking event and send. Newer filters such as suspension type track events in the format
 * {payloadKey: payloadValue} whereas older filters such as author track events in the format
 * {"payload": payloadValue}.
 * @param payloadKey Payload field, if any. Defaults to "payload" if not specified.
 * @param payloadValue Value to send as payload.
 * @returns
 */
export function buildAnalyticsPayload(
  payloadValue: string,
  payloadKey?: string
): Record<string, unknown> {
  return payloadKey
    ? { [payloadKey]: payloadValue }
    : { payload: payloadValue };
}
