import { MenuItem } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import {
  CATEGORY_KEY,
  OnFilterFn,
  SelectCategoryValueView,
} from "src/components/common/Filter/common/entities";
import { SelectionIcon } from "../../../../common/style";
import { MenuItemWrapper } from "./style";

interface Props {
  categoryKey: CATEGORY_KEY;
  isMultiselect: boolean;
  menuItems: SelectCategoryValueView[];
  onFilter: OnFilterFn;
}

export default function FilterMenuItems({
  categoryKey,
  isMultiselect,
  menuItems,
  onFilter,
}: Props): JSX.Element {
  return (
    <>
      {menuItems.map(({ key, count, label, selected }) => (
        <MenuItemWrapper key={key} isSelected={selected}>
          <MenuItem
            icon={
              <SelectionIcon
                icon={selected ? IconNames.TICK : IconNames.BLANK}
              />
            }
            labelElement={count}
            onClick={() => onFilter(categoryKey, key)}
            shouldDismissPopover={!isMultiselect}
            text={label}
          />
        </MenuItemWrapper>
      ))}
    </>
  );
}
