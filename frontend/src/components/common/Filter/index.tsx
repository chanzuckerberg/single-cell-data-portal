import { ChangeEvent, ReactNode } from "react";
import { CategoryFilterId } from "src/common/hooks/useCategoryFilter";
import {
  CategoryValueKey,
  CategoryView,
  OnFilterFn,
  OntologyCategoryView,
  OntologyMultiPanelCategoryView,
  Range,
  RangeCategoryView,
  SelectCategoryValueView,
  SelectCategoryView,
  SetSearchValueFn,
} from "src/components/common/Filter/common/entities";
import { formatNumberToScale } from "src/components/common/Filter/common/utils";
import BasicFilter from "src/components/common/Filter/components/BasicFilter";
import FilterLabel from "src/components/common/Filter/components/FilterLabel";
import FilterMenu from "src/components/common/Filter/components/FilterMenu";
import { MAX_DISPLAYABLE_MENU_ITEMS } from "src/components/common/Filter/components/FilterMenu/style";
import FilterRange from "src/components/common/Filter/components/FilterRange";
import { FilterSearchState } from "src/components/common/Filter/components/FilterSearch/common/useFilterSearch";
import FilterViews from "src/components/common/Filter/components/FilterViews";
import FilterMultiPanelCategoryView from "src/components/common/Filter/components/FilterViews/components/FilterMultiPanelCategoryView";
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
            Content={(filterSearchState: FilterSearchState): JSX.Element =>
              buildBasicFilterContent(categoryView, onFilter, filterSearchState)
            }
            flipEnabled={categoryViews.label !== "Tissue (Ontology)"} // TODO(cc) review use of flipEnabled prop
            isDisabled={isDisabled}
            key={categoryView.label}
            tags={<FilterTags tags={buildFilterTags(categoryView, onFilter)} />}
            target={buildFilterLabel(label, isDisabled, tooltip)}
          />
        );
      })}
    </>
  );
}

/**
 * Build content model of basic filter depending on category type.
 * @param categoryView - View model of category to display.
 * @param onFilter - Function to execute on select of category value or remove of selected category value.
 * @param filterSearchState - Filter search value and corresponding functions to clear or update search value.
 * @returns Element representing content to display inside basic filter menu.
 */
function buildBasicFilterContent(
  categoryView: CategoryView,
  onFilter: OnFilterFn,
  filterSearchState: FilterSearchState
): JSX.Element {
  const { key } = categoryView;
  const { clearSearchValueFn, searchValue, setSearchValue } = filterSearchState;

  // Update onFilter function with clear search value function.
  const onFilterWithClearSearch = (
    categoryKey: CategoryFilterId,
    key: CategoryValueKey | null, // null for ranges.
    values: CategoryValueKey[] | Range
  ) => {
    onFilter(categoryKey, key, values);
    clearSearchValueFn();
  };

  // Handle ontology categories.
  if (isOntologyCategoryView(categoryView)) {
    return (
      <FilterViews
        categoryKey={key}
        isSearchable={categoryView.isSearchable}
        isZerosVisible={categoryView.isZerosVisible}
        onFilter={onFilter}
        onUpdateSearchValue={onUpdateSearchValue}
        views={categoryView.views}
      />
    );
  }

  // Handle select categories
  if (isSelectCategoryView(categoryView)) {
    const { pinnedValues, unpinnedValues, values } = categoryView;
    return (
      <FilterMenu
        categoryKey={key}
        isMultiselect // Can possibly be single select with future filter types
        isSearchable={values.length > MAX_DISPLAYABLE_MENU_ITEMS}
        onFilter={onFilter}
        onUpdateSearchValue={onUpdateSearchValue}
        pinnedValues={filterCategoryValuesWithCount(pinnedValues)}
        searchValue={searchValue}
        setSearchValue={setSearchValue}
        unpinnedValues={filterCategoryValuesWithCount(unpinnedValues)}
        values={filterCategoryValuesWithCount(values)}
      />
    );
  }

  // Handle ontology multi panel categories
  if (isOntologyMultiPanelCategoryView(categoryView)) {
    return (
      <FilterMultiPanelCategoryView
        categoryView={categoryView}
        isSearchable={isFilterMultiPanelSearchable(categoryView)}
        onFilter={onFilterWithClearSearch}
        searchValue={searchValue}
        setSearchValue={setSearchValue}
      />
    );
  }

  // Otherwise, handle range categories
  return <FilterRange categoryView={categoryView} onFilter={onFilter} />;
}

/**
 * Build the filter label for the given category.
 * @param label - Category view label.
 * @param isDisabled - True if this category view is currently disabled.
 * @param tooltip - Category tooltip.
 * @returns React node representing content to display as filter label.
 */
function buildFilterLabel(
  label: string,
  isDisabled: boolean,
  tooltip?: string
): ReactNode {
  return (
    <FilterLabel isDisabled={isDisabled} label={label} tooltip={tooltip} />
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
 * @param categoryKey
 * @param onFilter
 * @returns ontology category tag.
 */
function buildOntologyCategoryTags(
  categoryView: OntologyCategoryView,
  categoryKey: CategoryFilterId,
  onFilter: OnFilterFn
): CategoryTag[] | undefined {
  return categoryView.views?.reduce((accum, species) => {
    species.selectedViews.forEach(({ key, label, values }) => {
      accum.push({
        label: label,
        onRemove: () => onFilter(categoryKey, key, values),
      });
    });
    return accum;
  }, [] as CategoryTag[]);
}

/**
 * Returns ontology multi panel category tag with tag label and corresponding Tag onRemove function.
 * TODO(cc) review code.
 * @param categoryView - Ontology multi panel category view.
 * @param categoryKey
 * @param onFilter - Function to execute on select of category value or remove of selected category value.
 * @returns ontology multi panel category tag.
 */
function buildOntologyMultiPanelCategoryTags(
  categoryView: OntologyMultiPanelCategoryView,
  categoryKey: CategoryFilterId,
  onFilter: OnFilterFn
): CategoryTag[] | undefined {
  const { panels } = categoryView;
  return panels.reduce((accum, ontologyCategoryView) => {
    ontologyCategoryView.views.forEach(({ key, label, selected, values }) => {
      if (selected) {
        accum.push({
          label: label,
          onRemove: () => onFilter(categoryKey, key, values),
        });
      }
    });
    return accum;
  }, [] as CategoryTag[]);
}

/**
 * Returns range category tag with tag label (the selected range) and corresponding Tag onRemove function.
 * @param categoryView
 * @param categoryKey
 * @param onFilter
 * @returns range category tag.
 */
function buildRangeCategoryTag(
  categoryView: RangeCategoryView,
  categoryKey: CategoryFilterId,
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
        // eslint-disable-next-line @typescript-eslint/ban-ts-comment -- -- TODO(cc) revisit clear - how to indicate nothing?
        // @ts-ignore
        onRemove: () => onFilter(categoryKey, null, {}),
      },
    ];
  }
}

/**
 * Returns selected category tags with tag label (the selected metadata label) and corresponding Tag onRemove function.
 * @param categoryView
 * @param categoryKey
 * @param onFilter
 * @returns selected category tags.
 */
function buildSelectCategoryTags(
  categoryView: SelectCategoryView,
  categoryKey: CategoryFilterId,
  onFilter: OnFilterFn
): CategoryTag[] {
  const { values } = categoryView;
  return values
    .filter((value) => value.selected)
    .map(({ key, label, values }) => {
      return {
        label: label,
        onRemove: () => onFilter(categoryKey, key, values),
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
 * Returns filtered category values where category count is greater than zero.
 * @param categoryValues - Category value view models for a given category.
 * @returns category values with a count
 */
function filterCategoryValuesWithCount(
  categoryValues: SelectCategoryValueView[]
): SelectCategoryValueView[] {
  return categoryValues.filter(({ count }) => count > 0);
}

/**
 * Returns true if ontology multi panel is searchable.
 * @param categoryView - Ontology multi panel category view.
 * @returns true if ontology multi panel is searchable.
 */
function isFilterMultiPanelSearchable(
  categoryView: OntologyMultiPanelCategoryView
): boolean {
  let isSearchable = false;
  for (const panel of categoryView.panels) {
    const { views } = panel;
    if (views.length > MAX_DISPLAYABLE_MENU_ITEMS) {
      isSearchable = true;
      break;
    }
  }
  return isSearchable;
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
 * Determine if the given category view is a range category view and not a select or ontology with or without multi panel category view.
 * TODO(cc) unused function.
 * @param categoryView - Selected filter value, either a category value key (e.g. "normal"), range (e.g. [0, 10]) or
 * ontology tree.
 * @returns True if the given category view is a range category view.
 */
export function isRangeCategoryView(
  categoryView: CategoryView
): categoryView is RangeCategoryView {
  return (categoryView as RangeCategoryView).max !== undefined;
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

/**
 * Sets state searchValue with updated search value.
 * @param changeEvent
 * @param setSearchValue
 */
function onUpdateSearchValue(
  changeEvent: ChangeEvent<HTMLInputElement>,
  setSearchValue: SetSearchValueFn
): void {
  setSearchValue(changeEvent.target.value.toLowerCase());
}
