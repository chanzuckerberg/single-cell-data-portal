import { InputGroup, Menu, MenuItem } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { Fragment, useEffect, useRef, useState } from "react";
import {
  CategoryValueView,
  CATEGORY_KEY,
  FilterCategoryValuesFn,
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
  multiselect: boolean;
  onFilter: OnFilterFn;
  onUpdateSearchValue: OnUpdateSearchValueFn;
  searchable: boolean;
  values: CategoryValueView[];
}

const ADDITIONAL_MENU_WIDTH = 8;

export default function FilterMenu({
  categoryKey,
  filterCategoryValues,
  multiselect,
  onFilter,
  onUpdateSearchValue,
  searchable,
  values,
}: Props): JSX.Element {
  const menuRef = useRef<HTMLSpanElement>(null);
  const [menuWidth, setMenuWidth] = useState(0);
  const [searchValue, setSearchValue] = useState<string>("");
  const menuItems = searchable
    ? filterCategoryValues(values, searchValue)
    : values;
  const emptyItems = menuItems.length === 0;
  const scrollable = menuItems.length > MAX_DISPLAYABLE_MENU_ITEMS;
  const MenuItemsScroller =
    searchable && scrollable ? MenuItemsWrapper : Fragment;

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
