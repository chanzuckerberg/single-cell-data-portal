import { MenuItem } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import {
  CATEGORY_FILTER_ID,
  OnFilterFn,
  SelectCategoryValueView,
} from "src/components/common/Filter/common/entities";
import { SelectionIcon } from "../../../../../../common/style";
import { MenuItemWrapper } from "./style";

interface Props {
  categoryFilterId: CATEGORY_FILTER_ID;
  isMultiselect: boolean;
  menuItems: SelectCategoryValueView[];
  onFilter: OnFilterFn;
}

export default function FilterMenuItems({
  categoryFilterId,
  isMultiselect,
  menuItems,
  onFilter,
}: Props): JSX.Element {
  return (
    <>
      {menuItems.map(({ key, count, label, selected, value }) => (
        <MenuItemWrapper key={key} isSelected={selected}>
          <MenuItem
            icon={
              <SelectionIcon
                icon={selected ? IconNames.TICK : IconNames.BLANK}
              />
            }
            labelElement={count}
            onClick={() => onFilter(categoryFilterId, value)}
            shouldDismissPopover={!isMultiselect}
            text={label}
          />
        </MenuItemWrapper>
      ))}
    </>
  );
}
