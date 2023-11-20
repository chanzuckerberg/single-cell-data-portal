import { Icon, ListItem } from "@czi-sds/components";
import React, { Fragment, ReactElement } from "react";
import { ListItemButton, ListItemIcon, ListItemText } from "@mui/material";
import {
  CATEGORY_FILTER_ID,
  OnFilterFn,
  OntologyCategoryTreeNodeView,
  SelectCategoryValueView,
} from "src/components/common/Filter/common/entities";
import {
  List,
  NoMatches,
} from "src/components/common/Filter/components/FilterContent/components/common/style";
import { ViewSublist } from "src/components/common/Filter/components/FilterContent/components/FilterViews/components/FilterViewList/style";

interface Props {
  categoryFilterId: CATEGORY_FILTER_ID;
  isZerosVisible: boolean;
  nested?: boolean;
  onFilter: OnFilterFn;
  values: OntologyCategoryTreeNodeView[] | SelectCategoryValueView[];
  ViewHeader?: ReactElement;
}

/**
 * Returns metadata values with a count, or all metadata values if zero count is specified as visible.
 * @param values - Metadata values for single or multiselect category.
 * @param isZerosVisible - Metadata value with a zero count are visible.
 * @returns metadata values for single or multiselect categories with a count and values without a count if specified.
 */
function filterCategoryValues(
  values: OntologyCategoryTreeNodeView[] | SelectCategoryValueView[],
  isZerosVisible: boolean
): OntologyCategoryTreeNodeView[] | SelectCategoryValueView[] {
  if (isZerosVisible) {
    return values;
  }
  return values.filter((value) => value.count);
}

/**
 * Determine if the given category value is an ontology category value and not a select value.
 * @param value - Metadata value for single or multiselect category.
 * @returns True if the given category value is an ontology category value.
 */
function isOntologyCategoryTreeNodeView(
  value: OntologyCategoryTreeNodeView | SelectCategoryValueView
): value is OntologyCategoryTreeNodeView {
  return (value as OntologyCategoryTreeNodeView).children !== undefined;
}

export default function FilterViewList({
  categoryFilterId,
  isZerosVisible,
  nested = false,
  onFilter,
  values,
  ViewHeader = undefined,
}: Props): JSX.Element {
  const filteredValues = filterCategoryValues(values, isZerosVisible);
  const ViewList = nested ? ViewSublist : List;
  return (
    <ViewList dense disablePadding subheader={ViewHeader}>
      {/* No matches */}
      {filteredValues.length === 0 ? (
        <NoMatches>No items found</NoMatches>
      ) : (
        filteredValues.map((filteredValue) => {
          const { categoryValueId, count, label, selected, selectedPartial } =
            filteredValue;
          let children;
          if (isOntologyCategoryTreeNodeView(filteredValue)) {
            children = filteredValue.children;
          }
          return (
            <Fragment key={categoryValueId}>
              {/* List item */}
              <ListItem>
                <ListItemButton
                  disabled={!count}
                  onClick={() =>
                    onFilter(categoryFilterId, categoryValueId, label)
                  }
                  selected={selected || selectedPartial}
                >
                  {/* Icon */}
                  <ListItemIcon>
                    {(selected || selectedPartial) && (
                      <Icon
                        sdsIcon={selected ? "check" : "minus"}
                        sdsSize="s"
                        sdsType="iconButton"
                      />
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
              {/* Nested list */}
              {children && !!children.length && (
                <FilterViewList
                  categoryFilterId={categoryFilterId}
                  isZerosVisible={isZerosVisible}
                  nested
                  onFilter={onFilter}
                  values={children}
                />
              )}
            </Fragment>
          );
        })
      )}
    </ViewList>
  );
}
