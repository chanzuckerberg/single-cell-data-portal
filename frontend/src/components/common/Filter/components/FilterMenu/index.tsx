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
  multiselect: boolean;
  onFilter: OnFilterFn;
  onUpdateSearchValue: OnUpdateSearchValueFn;
  searchable: boolean;
  values: CategoryValueView[];
}

export default function FilterMenu({
  categoryKey,
  filterCategoryValues,
  filterCategoryValuesWithCount,
  multiselect,
  onFilter,
  onUpdateSearchValue,
  searchable,
  values,
}: Props): JSX.Element {
  const menuRef = useRef<HTMLSpanElement>(null);
  const [menuWidth, setMenuWidth] = useState(0);
  const [searchValue, setSearchValue] = useState<string>("");
  const menuItemsWithCount = filterCategoryValuesWithCount(values);
  const menuItems =
    searchable && searchValue
      ? filterCategoryValues(menuItemsWithCount, searchValue)
      : menuItemsWithCount;
  const emptyItems = menuItems.length === 0;
  const scrollable = menuItems.length > MAX_DISPLAYABLE_MENU_ITEMS;
  const MenuItemsScroller =
    searchable && scrollable ? MenuItemsWrapper : Fragment;

  // Set initial width on menu to prevent resizing on filter of menu items.
  useEffect(() => {
    if (menuRef.current) {
      setMenuWidth(menuRef.current.children[0]?.clientWidth);
    }
  }, []);

  return (
    <MenuWrapper ref={menuRef} width={menuWidth}>
      <Menu>
        {searchable ? (
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
          <NoMatches shouldDismissPopover={false} text={"No matches found"} />
        ) : (
          <MenuItemsScroller>
            {menuItems.map(({ key, count, label, selected }) => (
              <MenuItemWrapper key={key} isSelected={selected}>
                <MenuItem
                  icon={selected ? IconNames.TICK : IconNames.BLANK}
                  labelElement={count}
                  onClick={() => onFilter(categoryKey, key)}
                  shouldDismissPopover={!multiselect}
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
