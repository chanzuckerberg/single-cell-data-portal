import { ChangeEvent } from "react";
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
import { MAX_DISPLAYABLE_MENU_ITEMS } from "src/components/common/Filter/components/FilterMenu/style";
import FilterMultiPanel from "src/components/common/Filter/components/FilterMultiPanel";
import FilterRange from "src/components/common/Filter/components/FilterRange";
import BasicFilter from "./components/BasicFilter";
import FilterLabel from "./components/FilterLabel";
import FilterMenu from "./components/FilterMenu";
import FilterTags, { CategoryTag } from "./components/FilterTags";

interface Props {
  categories: CategoryView[];
  onFilter: OnFilterFn;
}

export default function Filter({ categories, onFilter }: Props): JSX.Element {
  return (
    <>
      {categories.map((categoryView: CategoryView) => {
        const { key, label } = categoryView;
        const { values } = categoryView as SelectCategoryView;
        const { species } = categoryView as OntologyCategoryView;
        const isDisabled = isCategoryNA(categoryView); // TODO(cc) check isOntologyCategoryViewNA method return value
        return (
          <BasicFilter
            content={
              isOntologyCategoryView(categoryView) ? (
                <FilterMultiPanel
                  categoryKey={key}
                  onFilter={onFilter}
                  species={species}
                />
              ) : isSelectCategoryView(categoryView) ? (
                <FilterMenu
                  categoryKey={key}
                  filterCategoryValues={filterCategoryValues}
                  filterCategoryValuesWithCount={filterCategoryValuesWithCount}
                  isMultiselect // Can possibly be single select with future filter types
                  isSearchable={values.length > MAX_DISPLAYABLE_MENU_ITEMS}
                  onFilter={onFilter}
                  onUpdateSearchValue={onUpdateSearchValue}
                  values={values}
                />
              ) : (
                <FilterRange categoryView={categoryView} onFilter={onFilter} />
              )
            }
            isDisabled={isDisabled}
            key={key}
            tags={
              <FilterTags
                tags={
                  isSelectCategoryView(categoryView)
                    ? buildSelectCategoryTags(categoryView, key, onFilter)
                    : buildRangeCategoryTag(categoryView, key, onFilter)
                }
              />
            }
            target={<FilterLabel isDisabled={isDisabled} label={label} />}
          />
        );
      })}
    </>
  );
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
): CategoryTag[] | undefined {
  const { values } = categoryView;
  if (!values) {
    return;
  }
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
 * Returns filtered category values where category key includes search value.
 * @param values
 * @param searchValue
 * @returns array of category values filtered by the given search value
 */
function filterCategoryValues(
  values: SelectCategoryValueView[],
  searchValue: string
): SelectCategoryValueView[] {
  return values.filter(({ key }) => key.toLowerCase().includes(searchValue));
}

/**
 * Returns filtered category values where category count is greater than zero.
 * @param values
 * @returns category values with a count
 */
function filterCategoryValuesWithCount(
  values: SelectCategoryValueView[]
): SelectCategoryValueView[] {
  return values.filter(({ count }) => count > 0);
}

/**
 * Returns true if category is not applicable.
 * @param categoryView
 * @returns true when category  is not applicable.
 */
function isCategoryNA(categoryView: CategoryView): boolean {
  if (isSelectCategoryView(categoryView)) {
    return isSelectCategoryNA(categoryView);
  }
  if (isOntologyCategoryView(categoryView)) {
    return isOntologyCategoryViewNA(categoryView);
  }
  return isRangeCategoryNA(categoryView);
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
  return (categoryView as OntologyCategoryView).species !== undefined;
}

/**
 * Returns true if ontology category is not applicable, that is, there are species ontology trees.
 * @param categoryView
 * @returns true when range min and max are both 0 or both equal. // TODO(cc) review return statement
 */
function isOntologyCategoryViewNA(categoryView: OntologyCategoryView): boolean {
  return !categoryView.species || categoryView.species.length === 0; // TODO(cc) review return
}

/**
 * Returns true if range category is not applicable, that is, range min and max are both 0 or both equal.
 * @param categoryView
 * @returns true when range min and max are both 0 or both equal.
 */
function isRangeCategoryNA(categoryView: RangeCategoryView): boolean {
  const { max, min } = categoryView;
  return (min === 0 && max === 0) || min === max;
}

/**
 * Returns true if select category is not applicable, that is, all values have a count of 0.
 * @param categoryView
 * @returns true when all category values have a count of 0.
 */
function isSelectCategoryNA(categoryView: SelectCategoryView): boolean {
  const { values } = categoryView;
  return values?.every((value) => value.count === 0);
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
