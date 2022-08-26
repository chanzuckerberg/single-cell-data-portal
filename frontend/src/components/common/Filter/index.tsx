import {
  CategoryView,
  CATEGORY_FILTER_ID,
  OnFilterFn,
  OntologyCategoryView,
  OntologyMultiPanelCategoryView,
  RangeCategoryView,
  SelectCategoryView,
} from "src/components/common/Filter/common/entities";
import { formatNumberToScale } from "src/components/common/Filter/common/utils";
import BasicFilter from "src/components/common/Filter/components/BasicFilter";
import FilterLabel from "src/components/common/Filter/components/FilterLabel";
import FilterContent from "./components/FilterContent";
import FilterTags, { CategoryTag } from "./components/FilterTags";

interface Props {
  categoryViews: CategoryView[];
  onFilter: OnFilterFn;
}

export default function Filter({
  categoryViews,
  onFilter,
}: Props): JSX.Element {
  return (
    <>
      {categoryViews.map((categoryView: CategoryView) => {
        const { isDisabled = false, label, tooltip } = categoryView;
        return (
          <BasicFilter
            content={
              <FilterContent categoryView={categoryView} onFilter={onFilter} />
            }
            isDisabled={isDisabled}
            key={categoryView.label}
            tags={<FilterTags tags={buildFilterTags(categoryView, onFilter)} />}
            target={
              <FilterLabel
                isDisabled={isDisabled}
                label={label}
                tooltip={tooltip}
              />
            }
          />
        );
      })}
    </>
  );
}

/**
 * Build up the set of selected filter tags depending on category type.
 * @param categoryView - View model of category to display.
 * @param onFilter - Function to execute on select of category value or remove of selected category value.
 * @returns Array of category tags to be displayed as selected tags.
 */
function buildFilterTags(
  categoryView: CategoryView,
  onFilter: OnFilterFn
): CategoryTag[] | undefined {
  const { key } = categoryView;

  // Handle ontology categories
  if (isOntologyCategoryView(categoryView)) {
    return buildOntologyCategoryTags(categoryView, key, onFilter);
  }

  // Handle select categories
  if (isSelectCategoryView(categoryView)) {
    return buildSelectCategoryTags(categoryView, key, onFilter);
  }

  // Handle ontology multi panel categories
  if (isOntologyMultiPanelCategoryView(categoryView)) {
    return buildOntologyMultiPanelCategoryTags(categoryView, key, onFilter);
  }

  // Otherwise, handle range categories
  return buildRangeCategoryTag(categoryView, key, onFilter);
}

/**
 * Returns ontology category tag with tag label and corresponding Tag onRemove function.
 * @param categoryView
 * @param categoryFilterId
 * @param onFilter
 * @returns ontology category tag.
 */
function buildOntologyCategoryTags(
  categoryView: OntologyCategoryView,
  categoryFilterId: CATEGORY_FILTER_ID,
  onFilter: OnFilterFn
): CategoryTag[] | undefined {
  return categoryView.views?.reduce((accum, species) => {
    species.selectedViews.forEach(({ key, label, value }) => {
      accum.push({
        label: label,
        onRemove: () => onFilter(categoryFilterId, key, value),
      });
    });
    return accum;
  }, [] as CategoryTag[]);
}

/**
 * Returns ontology multi panel category tag with tag label and corresponding Tag onRemove function.
 * TODO(cc) review code.
 * @param categoryView - Ontology multi panel category view.
 * @param categoryFilterId
 * @param onFilter - Function to execute on select of category value or remove of selected category value.
 * @returns ontology multi panel category tag.
 */
function buildOntologyMultiPanelCategoryTags(
  categoryView: OntologyMultiPanelCategoryView,
  categoryFilterId: CATEGORY_FILTER_ID,
  onFilter: OnFilterFn
): CategoryTag[] | undefined {
  const { panels } = categoryView;
  return panels.reduce((accum, ontologyCategoryView) => {
    ontologyCategoryView.views.forEach(({ key, label, selected, value }) => {
      if (selected) {
        accum.push({
          label: label,
          onRemove: () => onFilter(categoryFilterId, key, value),
        });
      }
    });
    return accum;
  }, [] as CategoryTag[]);
}

/**
 * Returns range category tag with tag label (the selected range) and corresponding Tag onRemove function.
 * @param categoryView
 * @param categoryFilterId
 * @param onFilter
 * @returns range category tag.
 */
function buildRangeCategoryTag(
  categoryView: RangeCategoryView,
  categoryFilterId: CATEGORY_FILTER_ID,
  onFilter: OnFilterFn
): CategoryTag[] | undefined {
  const { selectedMax, selectedMin } = categoryView;
  if (!selectedMin && !selectedMax) {
    return;
  }
  if (selectedMin && selectedMax) {
    // There will only ever be a single selected tag for a range category but tag component is expecting an array:
    // create singleton array.
    return [
      {
        label: createRangeTagLabel(selectedMin, selectedMax),
        onRemove: () => onFilter(categoryFilterId, null, []),
      },
    ];
  }
}

/**
 * Returns selected category tags with tag label (the selected metadata label) and corresponding Tag onRemove function.
 * @param categoryView
 * @param categoryFilterId
 * @param onFilter
 * @returns selected category tags.
 */
function buildSelectCategoryTags(
  categoryView: SelectCategoryView,
  categoryFilterId: CATEGORY_FILTER_ID,
  onFilter: OnFilterFn
): CategoryTag[] {
  const { values } = categoryView;
  return values
    .filter((value) => value.selected)
    .map(({ key, label, value }) => {
      return {
        label: label,
        onRemove: () => onFilter(categoryFilterId, key, value),
      };
    });
}

/**
 * Returns filter tag label for the selected range of the slider.
 * @param min
 * @param max
 * @returns string portraying the selected range of the slider.
 */
function createRangeTagLabel(min: number, max: number): [string, string] {
  const minLabel = formatNumberToScale(min);
  const maxLabel = formatNumberToScale(max);
  return [minLabel, maxLabel];
}

/**
 * Determine if the given category view is an ontology category view and not a select or range or ontology multi panel category view.
 * @param categoryView - Selected filter value, either a category value key (e.g. "normal"), range (e.g. [0, 10]) or
 * ontology tree.
 * @returns True if the given category view is a select category view.
 */
export function isOntologyCategoryView(
  categoryView: CategoryView
): categoryView is OntologyCategoryView {
  return (categoryView as OntologyCategoryView).views !== undefined;
}

/**
 * Determine if the given category view is an ontology multi panel category view and not a select or range or ontology category view.
 * @param categoryView
 */
export function isOntologyMultiPanelCategoryView(
  categoryView: CategoryView
): categoryView is OntologyMultiPanelCategoryView {
  return (categoryView as OntologyMultiPanelCategoryView).panels !== undefined;
}

/**
 * Determine if the given category view is a selected category view and not an ontology with or without multi panel or range category view.
 * @param categoryView - Selected filter value, either a category value key (e.g. "normal"), range (e.g. [0, 10]) or
 * ontology tree.
 * @returns True if the given category view is a select category view.
 */
export function isSelectCategoryView(
  categoryView: CategoryView
): categoryView is SelectCategoryView {
  return (categoryView as SelectCategoryView).values !== undefined;
}
