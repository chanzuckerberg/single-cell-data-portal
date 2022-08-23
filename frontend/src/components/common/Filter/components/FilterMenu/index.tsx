import { Menu } from "@blueprintjs/core";
import { Fragment, useEffect, useRef, useState } from "react";
import { CategoryFilterId } from "src/common/hooks/useCategoryFilter";
import {
  OnFilterFn,
  OnUpdateSearchValueFn,
  SelectCategoryValueView,
  SetSearchValueFn,
} from "src/components/common/Filter/common/entities";
import FilterMenuItems from "src/components/common/Filter/components/FilterMenu/components/FilterMenuItems";
import FilterSearch from "src/components/common/Filter/components/FilterSearch";
import {
  MAX_DISPLAYABLE_MENU_ITEMS,
  MenuDivider,
  MenuItemsWrapper,
  MenuWrapper,
  NoMatches,
} from "./style";

interface Props {
  categoryKey: CategoryFilterId;
  isMultiselect: boolean;
  isSearchable: boolean;
  onFilter: OnFilterFn;
  onUpdateSearchValue: OnUpdateSearchValueFn;
  pinnedValues: SelectCategoryValueView[];
  searchValue: string;
  setSearchValue: SetSearchValueFn;
  unpinnedValues: SelectCategoryValueView[];
  values: SelectCategoryValueView[];
}

const ADDITIONAL_MENU_WIDTH = 8;

export default function FilterMenu({
  categoryKey,
  isMultiselect,
  isSearchable,
  onFilter,
  pinnedValues,
  searchValue,
  setSearchValue,
  unpinnedValues,
  values,
}: Props): JSX.Element {
  const menuRef = useRef<HTMLSpanElement>(null);
  const [menuWidth, setMenuWidth] = useState(0);
  const filteredPinnedValues = filterCategoryValues(pinnedValues, searchValue);
  const filteredUnpinnedValues = filterCategoryValues(
    unpinnedValues,
    searchValue
  );
  const filteredValues = filterCategoryValues(values, searchValue);
  const emptyItems = filteredValues.length === 0;
  const scrollable = isMenuScrollable(filteredValues, isSearchable);
  const MenuItemsScroller = scrollable ? MenuItemsWrapper : Fragment;
  const isMenuDivided =
    filteredPinnedValues.length > 0 && filteredUnpinnedValues.length > 0;
  const scrollerProps = scrollable ? { isMenuDivided } : {};

  // Set initial width on menu to prevent resizing on filter of menu items.
  useEffect(() => {
    if (menuRef.current) {
      setMenuWidth(
        menuRef.current?.children[0]?.clientWidth + ADDITIONAL_MENU_WIDTH
      );
    }
  }, []);

  return (
    <MenuWrapper menuWidth={menuWidth} ref={menuRef}>
      <Menu>
        {/* Optional search bar */}
        {isSearchable && (
          <FilterSearch
            searchValue={searchValue}
            setSearchValue={setSearchValue}
          />
        )}
        {/* No matches */}
        {emptyItems ? (
          <NoMatches
            multiline /* required when no matches text length is longer than longest category value */
            shouldDismissPopover={false}
            text={"No matches found"}
          />
        ) : (
          <MenuItemsScroller {...scrollerProps}>
            {/* Pinned values */}
            {filteredPinnedValues.length > 0 && (
              <FilterMenuItems
                categoryKey={categoryKey}
                isMultiselect={isMultiselect}
                menuItems={filteredPinnedValues}
                onFilter={onFilter}
              />
            )}
            {/* Menu divider */}
            {isMenuDivided && <MenuDivider />}
            {/* Unpinned values */}
            <FilterMenuItems
              categoryKey={categoryKey}
              isMultiselect={isMultiselect}
              menuItems={filteredUnpinnedValues}
              onFilter={onFilter}
            />
          </MenuItemsScroller>
        )}
      </Menu>
    </MenuWrapper>
  );
}

/**
 * Returns filtered category values where category key includes search value.
 * @param categoryValues - Category value view models for a given category.
 * @param searchValue - Search string to filters category values.
 * @returns array of category values filtered by the given search value
 */
function filterCategoryValues(
  categoryValues: SelectCategoryValueView[],
  searchValue: string
): SelectCategoryValueView[] {
  if (searchValue) {
    return categoryValues.filter(({ key }) =>
      key.toLowerCase().includes(searchValue)
    );
  }
  return categoryValues;
}

/**
 * Returns true if the menu is searchable and exceeds the maximum number of displayable menu items.
 * @param categoryValues - Category value view models for a given category.
 * @param isSearchable - Category is searchable.
 */
function isMenuScrollable(
  categoryValues: SelectCategoryValueView[],
  isSearchable: boolean
): boolean {
  return isSearchable && categoryValues.length > MAX_DISPLAYABLE_MENU_ITEMS;
}
