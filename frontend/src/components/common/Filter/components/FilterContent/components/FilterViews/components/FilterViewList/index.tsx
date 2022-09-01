import { IconNames } from "@blueprintjs/icons";
import { List } from "czifui";
import React, { Fragment, ReactElement } from "react";
import {
  CATEGORY_FILTER_ID,
  OnFilterFn,
  OntologyCategoryTreeNodeView,
  SelectCategoryValueView,
} from "src/components/common/Filter/common/entities";
import { SelectionIcon } from "src/components/common/Filter/common/style";
import {
  NoMatches,
  ViewListItem,
  ViewListItemText,
  ViewSublist,
} from "src/components/common/Filter/components/FilterContent/components/FilterViews/components/FilterViewList/style";

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
  return (value as OntologyCategoryTreeNodeView).selectedPartial !== undefined;
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
          const { key, count, label, selected, selectedPartial, value } =
            filteredValue;
          let children;
          // TODO(cc) revisit type predicate
          if (isOntologyCategoryTreeNodeView(filteredValue)) {
            children = filteredValue.children;
          }
          return (
            <Fragment key={key}>
              {/* List item */}
              <ViewListItem
                button
                disabled={!count}
                onClick={() => onFilter(categoryFilterId, key, value)}
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
