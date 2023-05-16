import { Icon, ListItem } from "czifui";
import React from "react";
import { ListItemButton, ListItemIcon, ListItemText } from "@mui/material";
import {
  CATEGORY_FILTER_ID,
  OnFilterFn,
  SelectCategoryValueView,
} from "src/components/common/Filter/common/entities";

interface Props {
  categoryFilterId: CATEGORY_FILTER_ID;
  menuItems: SelectCategoryValueView[];
  onFilter: OnFilterFn;
}

export default function FilterMenuItems({
  categoryFilterId,
  menuItems,
  onFilter,
}: Props): JSX.Element {
  return (
    <>
      {menuItems.map(({ categoryValueId, count, label, selected }) => (
        <ListItem key={categoryValueId}>
          <ListItemButton
            onClick={() => onFilter(categoryFilterId, categoryValueId, label)}
            selected={selected}
          >
            {/* Icon */}
            <ListItemIcon>
              {selected && (
                <Icon sdsIcon="check" sdsSize="s" sdsType="iconButton" />
              )}
            </ListItemIcon>
            {/* List item text and count */}
            <ListItemText
              disableTypography
              primary={<span>{label}</span>}
              secondary={<span>{count}</span>}
            />
          </ListItemButton>
        </ListItem>
      ))}
    </>
  );
}
