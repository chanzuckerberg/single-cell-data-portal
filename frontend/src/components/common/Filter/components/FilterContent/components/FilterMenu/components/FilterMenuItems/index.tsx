import { Icon } from "@czi-sds/components";
import React from "react";
import { ListItemIcon, ListItemText } from "@mui/material";
import {
  CATEGORY_FILTER_ID,
  OnFilterFn,
  SelectCategoryValueView,
} from "src/components/common/Filter/common/entities";
import { ListItem } from "./style";

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
  // (clevercanary): For optimal rendering performance for the PUBLICATION_AUTHORS category filter, the SDS ListItem
  // and Mui ListItemButton components are substituted a standard li element.
  return (
    <>
      {menuItems.map(({ categoryValueId, count, label, selected }) => (
        <ListItem
          key={categoryValueId}
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
        </ListItem>
      ))}
    </>
  );
}
