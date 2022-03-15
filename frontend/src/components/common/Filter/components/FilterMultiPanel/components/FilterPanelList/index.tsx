import { IconNames } from "@blueprintjs/icons";
import { Box, List, ListItem, ListItemText } from "@material-ui/core";
import { Fragment, ReactElement } from "react";
import {
  CATEGORY_KEY,
  OnFilterFn,
  OntologyCategoryValueView,
} from "src/components/common/Filter/common/entities";
import { SelectionIcon } from "src/components/common/Filter/common/style";
import { useFilterPanelListStyles } from "src/components/common/Filter/components/FilterMultiPanel/components/FilterPanelList/style";
import { GRAY } from "src/components/common/theme";

interface Props {
  categoryKey: CATEGORY_KEY;
  nested?: boolean;
  onFilter: OnFilterFn;
  PanelHeader?: ReactElement;
  values: OntologyCategoryValueView[];
}

export default function FilterPanelList({
  categoryKey,
  nested = false,
  onFilter,
  PanelHeader = undefined,
  values,
}: Props): JSX.Element {
  const classes = useFilterPanelListStyles();
  return (
    <List
      classes={{ root: nested ? classes.sublist : undefined }}
      dense
      disablePadding
      subheader={PanelHeader}
    >
      {values.map(
        ({ key, children, count, label, selectedPartial, selected }) => (
          <Fragment key={key}>
            {/* List item */}
            <ListItem
              alignItems="flex-start"
              button
              classes={{ root: classes.listItem }}
              component="li"
              dense
              disabled={!count}
              disableGutters
              disableRipple
              onClick={() => onFilter(categoryKey, key)}
            >
              {/* Icon - bp icon to uphold ui consistency between filter menu and filter multi panel */}
              <SelectionIcon
                icon={
                  selected
                    ? IconNames.TICK
                    : selectedPartial
                    ? IconNames.MINUS
                    : IconNames.BLANK
                }
              />
              {/* List item text and count */}
              <ListItemText
                classes={{ root: classes.listItemText }}
                disableTypography
                primary={
                  <Box
                    component="span"
                    flex={1}
                    fontWeight={selected ? 500 : undefined}
                    mr={4}
                  >
                    {label}
                  </Box>
                }
                secondary={
                  <Box color={GRAY.A} component="span">
                    {count}
                  </Box>
                }
              />
            </ListItem>
            {/* Nested list */}
            {children && children.length && (
              <FilterPanelList
                categoryKey={categoryKey}
                nested
                onFilter={onFilter}
                values={children}
              />
            )}
          </Fragment>
        )
      )}
    </List>
  );
}
