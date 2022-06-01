import { ChangeEvent, ReactNode } from "react";
import { CategoryKey } from "src/common/hooks/useCategoryFilter";
import {
  CategoryView,
  OnFilterFn,
  OntologyCategoryView,
  RangeCategoryView,
  SelectCategoryValueView,
  SelectCategoryView,
  SetSearchValueFn,
} from "src/components/common/Filter/common/entities";
import { formatNumberToScale } from "src/components/common/Filter/common/utils";
import FilterLabel from "src/components/common/Filter/components/FilterLabel";
import FilterMenu from "src/components/common/Filter/components/FilterMenu";
import { MAX_DISPLAYABLE_MENU_ITEMS } from "src/components/common/Filter/components/FilterMenu/style";
import FilterMultiPanel from "src/components/common/Filter/components/FilterMultiPanel";
import FilterRange from "src/components/common/Filter/components/FilterRange";
import BasicFilter from "./components/BasicFilter";
import FilterTags, { CategoryTag } from "./components/FilterTags";

interface Props {
  categories: CategoryView[];
  onFilter: OnFilterFn;
}

export default function Filter({ categories, onFilter }: Props): JSX.Element {
  return (
    <>
      {categories.map((categoryView: CategoryView) => {
        const { isDisabled = false, key } = categoryView;
        return (
          <BasicFilter
            content={buildBasicFilterContent(categoryView, onFilter)}
            isDisabled={isDisabled}
            key={key}
            tags={<FilterTags tags={buildFilterTags(categoryView, onFilter)} />}
            target={buildFilterLabel(categoryView, isDisabled)}
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
 * @returns React node representing content to display inside basic filter menu.
 */
function buildBasicFilterContent(
  categoryView: CategoryView,
  onFilter: OnFilterFn
): ReactNode {
  const { key } = categoryView;

  // Handle ontology categories.
  if (isOntologyCategoryView(categoryView)) {
    return (
      <FilterMultiPanel
        categoryKey={key}
        onFilter={onFilter}
        species={categoryView.views}
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
        unpinnedValues={filterCategoryValuesWithCount(unpinnedValues)}
        values={filterCategoryValuesWithCount(values)}
      />
    );
  }

  // Otherwise, handle range categories
  return <FilterRange categoryView={categoryView} onFilter={onFilter} />;
}

/**
 * Build the filter label for the given category.
 * @param categoryView - View model of category to display.
 * @param isDisabled - True if this category is currently disabled.
 * @returns React node representing content to display as filter label.
 */
function buildFilterLabel(
  categoryView: CategoryView,
  isDisabled: boolean
): ReactNode {
  const { label, tooltip } = categoryView;

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
  categoryKey: CategoryKey,
  onFilter: OnFilterFn
): CategoryTag[] | undefined {
  return categoryView.views?.reduce((accum, species) => {
    species.selectedViews.forEach(({ key, label }) => {
      accum.push({ label: label, onRemove: () => onFilter(categoryKey, key) });
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
  categoryKey: CategoryKey,
  onFilter: OnFilterFn
): CategoryTag[] | undefined {
  const { selectedMax, selectedMin } = categoryView;
  if (!selectedMin && !selectedMax) {
    return;
  }
  if (selectedMin && selectedMax) {
    // There will only ever be a single selected tag for a range category but tag component is expecting an array: create
    // singleton array.
    return [
      {
        label: createRangeTagLabel(selectedMin, selectedMax),
        onRemove: () => onFilter(categoryKey, []),
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
  categoryKey: CategoryKey,
  onFilter: OnFilterFn
): CategoryTag[] {
  const { values } = categoryView;
  return values
    .filter((value) => value.selected)
    .map(({ key, label }) => {
      return { label: label, onRemove: () => onFilter(categoryKey, key) };
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
 * Determine if the given category view is an ontology category view and not a select or range category view.
 * @param categoryView - Selected filter value, either a category value key (e.g. "normal"), range (e.g. [0, 10]) or
 * ontology tree.
 * @returns True if the given category view is a select category view.
 */
function isOntologyCategoryView(
  categoryView: CategoryView
): categoryView is OntologyCategoryView {
  return (categoryView as OntologyCategoryView).views !== undefined;
}

/**
 * Determine if the given category view is a selected category view and not an ontology or range category view.
 * @param categoryView - Selected filter value, either a category value key (e.g. "normal"), range (e.g. [0, 10]) or
 * ontology tree.
 * @returns True if the given category view is a select category view.
 */
function isSelectCategoryView(
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
