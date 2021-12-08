import { ChangeEvent } from "react";
import {
  CategoryValueView,
  CategoryView,
  OnFilterFn,
  SetSearchValueFn,
} from "src/components/common/Filter/common/entities";
import { MAX_DISPLAYABLE_MENU_ITEMS } from "src/components/common/Filter/components/FilterMenu/style";
import BasicFilter from "./components/BasicFilter";
import FilterLabel from "./components/FilterLabel";
import FilterMenu from "./components/FilterMenu";
import FilterTags from "./components/FilterTags";

interface Props {
  categories: CategoryView[];
  onFilter: OnFilterFn;
}

export default function Filter({ categories, onFilter }: Props): JSX.Element {
  return (
    <>
      {categories.map(({ label, key, values }) => {
        const isDisabled = isCategoryNA(values);
        return (
          <BasicFilter
            content={
              <FilterMenu
                categoryKey={key}
                filterCategoryValues={filterCategoryValues}
                filterCategoryValuesWithCount={filterCategoryValuesWithCount}
                isMultiselect={true} // Can possibly be single select with future filter types
                isSearchable={values.length > MAX_DISPLAYABLE_MENU_ITEMS}
                onFilter={onFilter}
                onUpdateSearchValue={onUpdateSearchValue}
                values={values}
              />
            }
            isDisabled={isDisabled}
            key={key}
            tags={
              <FilterTags
                categoryKey={key}
                onFilter={onFilter}
                selectedValues={filterSelectedValues(values)}
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
 * Returns filtered category values where category key includes search value.
 * @param values
 * @param searchValue
 * @returns array of category values filtered by the given search value
 */
function filterCategoryValues(
  values: CategoryValueView[],
  searchValue: string
): CategoryValueView[] {
  return values.filter(({ key }) => key.toLowerCase().includes(searchValue));
}

/**
 * Returns filtered category values where category count is greater than zero.
 * @param values
 * @returns category values with a count
 */
function filterCategoryValuesWithCount(
  values: CategoryValueView[]
): CategoryValueView[] {
  return values.filter(({ count }) => count > 0);
}

/**
 * Returns selected category values.
 * @param values
 * @returns selected category values
 */
function filterSelectedValues(
  values: CategoryValueView[]
): CategoryValueView[] {
  return values.filter((value) => value.selected);
}

/**
 * Returns true if category is not applicable, that is, all values have a count of 0.
 * @param values
 * @returns true when all category values have a count of 0
 */
function isCategoryNA(values: CategoryValueView[]) {
  return values.every((value) => value.count === 0);
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
