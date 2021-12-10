import { InputGroup, Menu, MenuItem } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { Fragment, useEffect, useRef, useState } from "react";
import {
  CategoryValueView,
  CATEGORY_KEY,
  FilterCategoryValuesFn,
  FilterCategoryValuesWithCountFn,
  OnFilterFn,
  OnUpdateSearchValueFn,
} from "src/components/common/Filter/common/entities";
import {
  InputGroupWrapper,
  MAX_DISPLAYABLE_MENU_ITEMS,
  MenuItemsWrapper,
  MenuItemWrapper,
  MenuWrapper,
  NoMatches,
} from "./style";

interface Props {
  categoryKey: CATEGORY_KEY;
  filterCategoryValues: FilterCategoryValuesFn;
  filterCategoryValuesWithCount: FilterCategoryValuesWithCountFn;
  isMultiselect: boolean;
  isSearchable: boolean;
  onFilter: OnFilterFn;
  onUpdateSearchValue: OnUpdateSearchValueFn;
  values: CategoryValueView[];
}

const ADDITIONAL_MENU_WIDTH = 8;

export default function FilterMenu({
  categoryKey,
  filterCategoryValues,
  filterCategoryValuesWithCount,
  isMultiselect,
  onFilter,
  onUpdateSearchValue,
  isSearchable,
  values,
}: Props): JSX.Element {
  const menuRef = useRef<HTMLSpanElement>(null);
  const [menuWidth, setMenuWidth] = useState(0);
  const [searchValue, setSearchValue] = useState<string>("");
  const menuItemsWithCount = filterCategoryValuesWithCount(values);
  const menuItems =
    isSearchable && searchValue
      ? filterCategoryValues(menuItemsWithCount, searchValue)
      : menuItemsWithCount;
  const emptyItems = menuItems.length === 0;
  const scrollable = menuItems.length > MAX_DISPLAYABLE_MENU_ITEMS;
  const MenuItemsScroller =
    isSearchable && scrollable ? MenuItemsWrapper : Fragment;

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
        {isSearchable ? (
          <InputGroupWrapper>
            <InputGroup
              leftIcon={IconNames.SEARCH}
              onChange={(changeEvent) =>
                onUpdateSearchValue(changeEvent, setSearchValue)
              }
              placeholder="Search"
            />
          </InputGroupWrapper>
        ) : null}
        {emptyItems ? (
          <NoMatches
            multiline /* required when no matches text length is longer than longest category value */
            shouldDismissPopover={false}
            text={"No matches found"}
          />
        ) : (
          <MenuItemsScroller>
            {menuItems.map(({ key, count, label, selected }) => (
              <MenuItemWrapper key={key} isSelected={selected}>
                <MenuItem
                  icon={selected ? IconNames.TICK : IconNames.BLANK}
                  labelElement={count}
                  onClick={() => onFilter(categoryKey, key)}
                  shouldDismissPopover={!isMultiselect}
                  text={label}
                />
              </MenuItemWrapper>
            ))}
          </MenuItemsScroller>
        )}
      </Menu>
    </MenuWrapper>
  );
}
