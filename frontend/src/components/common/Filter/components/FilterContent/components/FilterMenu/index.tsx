import React, { useEffect, useRef, useState } from "react";
import { Divider } from "@mui/material";
import { List } from "@czi-sds/components";
import {
  CATEGORY_FILTER_ID,
  OnFilterFn,
  SelectCategoryValueView,
} from "src/components/common/Filter/common/entities";
import FilterMenuItems from "src/components/common/Filter/components/FilterContent/components/FilterMenu/components/FilterMenuItems";
import FilterSearch from "src/components/common/Filter/components/FilterSearch";
import { SetSearchValueFn } from "src/components/common/Filter/components/FilterSearch/common/useFilterSearch";
import { NoMatches } from "src/components/common/Filter/components/FilterContent/components/common/style";
import { FilterMenu as Menu, MAX_DISPLAYABLE_MENU_ITEMS } from "./style";

interface Props {
  categoryFilterId: CATEGORY_FILTER_ID;
  isSearchable: boolean;
  onFilter: OnFilterFn;
  pinnedValues: SelectCategoryValueView[];
  searchValue: string;
  setSearchValue: SetSearchValueFn;
  unpinnedValues: SelectCategoryValueView[];
  values: SelectCategoryValueView[];
}

// Additional menu width facilitates the rendering of selected menu items (where font-weight is "semibold")
// without any wrapping.
export const ADDITIONAL_MENU_WIDTH = 24;

export default function FilterMenu({
  categoryFilterId,
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
  const isMenuDivided =
    filteredPinnedValues.length > 0 && filteredUnpinnedValues.length > 0;

  // Set initial width on menu to prevent resizing on filter of menu items.
  useEffect(() => {
    if (menuRef.current) {
      setMenuWidth(
        menuRef.current?.children[0]?.clientWidth + ADDITIONAL_MENU_WIDTH
      );
    }
  }, []);

  return (
    <Menu
      isMenuDivided={isMenuDivided}
      menuWidth={menuWidth}
      ref={menuRef}
      scrollable={scrollable}
    >
      {/* Optional search bar */}
      {isSearchable && (
        <FilterSearch
          searchValue={searchValue}
          setSearchValue={setSearchValue}
        />
      )}
      <List>
        {/* No matches */}
        {emptyItems ? (
          <NoMatches>No matches found</NoMatches>
        ) : (
          <>
            {/* Pinned values */}
            {filteredPinnedValues.length > 0 && (
              <FilterMenuItems
                categoryFilterId={categoryFilterId}
                menuItems={filteredPinnedValues}
                onFilter={onFilter}
              />
            )}
            {/* Menu divider */}
            {isMenuDivided && <Divider />}
            {/* Unpinned values */}
            <FilterMenuItems
              categoryFilterId={categoryFilterId}
              menuItems={filteredUnpinnedValues}
              onFilter={onFilter}
            />
          </>
        )}
      </List>
    </Menu>
  );
}

/**
 * Returns filtered category values where category key includes search value.
 * @param categoryValues - Category value view models for a given category.
 * @param searchValue - Search string to filters category values.
 * @returns array of category values filtered by the given search value.
 */
function filterCategoryValues(
  categoryValues: SelectCategoryValueView[],
  searchValue: string
): SelectCategoryValueView[] {
  if (searchValue) {
    return categoryValues.filter(({ categoryValueId }) =>
      categoryValueId.toLowerCase().includes(searchValue.toLowerCase())
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
