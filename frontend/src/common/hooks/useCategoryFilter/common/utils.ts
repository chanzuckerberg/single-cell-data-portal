import { Filters } from "react-table";
import { removeOntologyTermIdPrefix } from "src/common/hooks/useCategoryFilter/common/multiPanelOntologyUtils";
import { COLLATOR_CASE_INSENSITIVE } from "src/components/common/Filter/common/constants";
import {
  Categories,
  CATEGORY_FILTER_ID,
  CategoryFilter,
  CategoryFilterConfig,
  CategorySetValue,
  CategoryValueId,
  PUBLICATION_DATE_LABELS,
  SelectCategoryValueView,
  SELF_REPORTED_ETHNICITY_UNSPECIFIED_LABEL,
} from "src/components/common/Filter/common/entities";

/**
 * General utils shared across the different category filter types.
 */

/**
 * Find and return the selected values for the given category.
 * @param categoryFilterId - ID of category to find selected filters of.
 * @param filters - Current set of selected category values (values) keyed by category (id).
 * @returns Array of filters
 */
export function getCategoryFilter<T extends Categories>(
  categoryFilterId: CATEGORY_FILTER_ID,
  filters: Filters<T>
): CategoryFilter | undefined {
  return filters.find((filter) => filter.id === categoryFilterId);
}

/**
 * Build the display value for the given category and category value. For ontology terms, look up corresponding labels.
 * @param config - Config model of category to build category value views for.
 * @param categoryValueId - Category value to display (e.g. "normal").
 * @param ontologyTermLabelsById - Set of ontology term labels keyed by term ID, used to determine labels for ontology
 * @returns String to display as a label for the given category and category value.
 */
export function buildCategoryValueLabel(
  config: CategoryFilterConfig,
  categoryValueId: CategoryValueId,
  ontologyTermLabelsById: Map<string, string>
): string {
  const { categoryFilterId } = config;
  if (categoryFilterId === CATEGORY_FILTER_ID.PUBLICATION_DATE_VALUES) {
    return PUBLICATION_DATE_LABELS[
      `LABEL_${categoryValueId}` as keyof typeof PUBLICATION_DATE_LABELS
    ];
  }

  if (
    categoryFilterId === CATEGORY_FILTER_ID.SELF_REPORTED_ETHNICITY &&
    !isSelfReportedEthnicitySpecified(categoryValueId)
  ) {
    return SELF_REPORTED_ETHNICITY_UNSPECIFIED_LABEL[
      categoryValueId as keyof typeof SELF_REPORTED_ETHNICITY_UNSPECIFIED_LABEL
    ];
  }

  // Look up labels for ontology term IDs.
  if (config.labelKind === "LOOKUP_LABEL_BY_TERM_ID") {
    let processedCategoryValueKey = categoryValueId;
    if (
      config.categoryFilterId === CATEGORY_FILTER_ID.TISSUE_CALCULATED ||
      config.categoryFilterId === CATEGORY_FILTER_ID.CELL_TYPE_CALCULATED
    ) {
      processedCategoryValueKey = removeOntologyTermIdPrefix(categoryValueId);
    }

    return (
      ontologyTermLabelsById.get(processedCategoryValueKey) ?? categoryValueId
    );
  }

  // Return all other category values as is.
  return categoryValueId;
}

/**
 * Determine if the given category value is a select category value (and not a range category value).
 * @param categorySetValue - Range or category set value.
 * @returns True if category set value is a set of category value keys.
 */
export function isCategoryValueIdSet(
  categorySetValue: CategorySetValue
): categorySetValue is Set<CategoryValueId> {
  return categorySetValue instanceof Set;
}

/**
 * Determine if the given self-reported ethnicity is considered unspecified (that is, na or unknown).
 * @param categoryValueKey - Self-Reported Ethnicity value to check if it's specified.
 * @returns True if self-reported ethnicity is either na or unknown.
 */
function isSelfReportedEthnicitySpecified(categoryValueKey: CategoryValueId) {
  const label =
    SELF_REPORTED_ETHNICITY_UNSPECIFIED_LABEL[
      categoryValueKey as keyof typeof SELF_REPORTED_ETHNICITY_UNSPECIFIED_LABEL
    ];
  return !label;
}

/**
 * Sort category value views by the given key, ascending. Sort can be by either the key (required for values such
 * as publication date) or label (required for values such as tissue that are backed by an ontology term ID).
 * @param key - Value to sort category value views by.
 * @returns Function that returns a number indicating sort precedence of cv0 vs cv1.
 */
export function sortCategoryValueViews(key: "categoryValueId" | "label") {
  return (
    cvv0: SelectCategoryValueView,
    cvv1: SelectCategoryValueView
  ): number => {
    return COLLATOR_CASE_INSENSITIVE.compare(cvv0[key], cvv1[key]);
  };
}
