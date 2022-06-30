import { IconNames } from "@blueprintjs/icons";
import { List } from "czifui";
import React, { Fragment, ReactElement } from "react";
import {
  CATEGORY_KEY,
  OnFilterFn,
  OntologyCategoryTreeNodeView,
} from "src/components/common/Filter/common/entities";
import { SelectionIcon } from "src/components/common/Filter/common/style";
import {
  NoMatches,
  ViewListItem,
  ViewListItemText,
  ViewSublist,
} from "src/components/common/Filter/components/FilterViews/components/FilterViewList/style";

interface Props {
  categoryKey: CATEGORY_KEY;
  isZerosVisible: boolean;
  nested?: boolean;
  onFilter: OnFilterFn;
  values: OntologyCategoryTreeNodeView[];
  ViewHeader?: ReactElement;
}

export default function FilterViewList({
  categoryKey,
  isZerosVisible,
  nested = false,
  onFilter,
  values,
  ViewHeader = undefined,
}: Props): JSX.Element {
  const filteredValues = values.filter(
    (value) => isZerosVisible || value.count
  );
  const ViewList = nested ? ViewSublist : List;
  return (
    <ViewList dense disablePadding subheader={ViewHeader}>
      {/* No matches */}
      {filteredValues.length === 0 ? (
        <NoMatches>No items found</NoMatches>
      ) : (
        filteredValues.map(
          ({ key, children, count, label, selectedPartial, selected }) => (
            <Fragment key={key}>
              {/* List item */}
              <ViewListItem
                button
                disabled={!count}
                onClick={() => onFilter(categoryKey, key)}
              >
                {/* Icon - bp icon to uphold ui consistency between filter menu and filter views */}
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
                <ViewListItemText
                  disableTypography
                  primary={<span>{label}</span>}
                  secondary={<span>{count}</span>}
                  selected={selected}
                />
              </ViewListItem>
              {/* Nested list */}
              {children && children.length && (
                <FilterViewList
                  categoryKey={categoryKey}
                  isZerosVisible={isZerosVisible}
                  nested
                  onFilter={onFilter}
                  values={children}
                />
              )}
            </Fragment>
          )
        )
      )}
    </ViewList>
  );
}
