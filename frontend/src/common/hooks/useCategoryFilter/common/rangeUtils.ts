import { track } from "src/common/analytics";
import { TOOLTIP_CATEGORY_DISABLED } from "src/components/common/Filter/common/constants";
import {
  CategoryFilterConfig,
  CategorySet,
  CategoryValueId,
  CATEGORY_FILTER_ID,
  FilterState,
  Range,
  RangeCategory,
  RangeCategoryView,
} from "src/components/common/Filter/common/entities";
import { isCategoryValueIdSet } from "./utils";

/**
 * Utils specific to range category filters.
 */

/**
 * Add range categories to the next filter state. Range categories are not summarized and must be explicitly added from
 * the original category set state.
 * @param categoryFilterIds - Set of category IDs to include for this filter instance.
 * @param nextFilterState - Filter state currently being built due to change in filter.
 * @param categorySet - Original, unfiltered sets of category values keyed by their category.
 */
export function addRangeCategories(
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
 * Build view model of range category.
 * @param config - Config model of ontology category.
 * @param rangeCategory - Internal filter model of range category.
 * @returns Range view model.
 */
export function buildRangeCategoryView(
  config: CategoryFilterConfig,
  rangeCategory: RangeCategory
): RangeCategoryView {
  const { categoryFilterId } = config;

  // Build view model of range category.
  const rangeView: RangeCategoryView = {
    categoryFilterId: categoryFilterId,
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
 * Returns true if range category is disabled, that is, range min and max are both 0 or both equal.
 * @param categoryView - Range category view to check enabled/disabled state of.
 * @returns true when range min and max are both 0 or both equal.
 */
function isRangeCategoryDisabled(categoryView: RangeCategoryView): boolean {
  const { max, min } = categoryView;
  return (min === 0 && max === 0) || min === max;
}

/**
 * Determine if the given selected value is a range and a category value ID.
 * @param selectedValue - Selected filter value, either a category value ID (e.g. "normal"), or a
 * range (e.g. [1,100]).
 * @returns True if given selected value is a range.
 */
export function isSelectedValueRange(
  selectedValue: CategoryValueId | Range
): selectedValue is Range {
  return Array.isArray(selectedValue);
}

/**
 * Handle select of range min/max value: set next set of filters for this category. Track updated range.
 * @param config - Configuration model of selected category.
 * @param selectedValue - Selected category value key (e.g. [1, 100]).
 * @returns The selected value with the latest range selection applied.
 */
export function onFilterRangeCategory(
  config: CategoryFilterConfig,
  selectedValue: Range
) {
  const { analyticsEvent } = config;

  // Track select of new range mim/max, ignoring any clear of selected range. Only track if event is specified on
  // configuration model
  if (analyticsEvent && selectedValue.length > 0) {
    const [min, max] = selectedValue;
    track(analyticsEvent, {
      max,
      min,
    });
  }

  // Return the filters for this range category.
  return selectedValue;
}
