import React, { Fragment } from "react";
import {
  CategoryView,
  OnFilterFn,
} from "src/components/common/Filter/common/entities";
import { MAX_DISPLAYABLE_LIST_ITEMS } from "src/components/common/Filter/components/FilterViews/components/FilterView";
import {
  ViewDivider,
  VIEW_LIST_ITEM_HEIGHT,
  VIEW_LIST_SUBHEADER_HEIGHT,
} from "src/components/common/Filter/components/FilterViews/components/FilterView/style";
import { ViewsMenu } from "src/components/common/Filter/components/FilterViews/style";
import FilterCategoryView from "../FilterCategoryView";

interface Props {
  categoryViews: CategoryView[];
  onFilter: OnFilterFn;
}

/**
 * Returns the maximum possible pixel height of any view list.
 * @param views - Array of category views objects.
 * @return maximum possible view list height in pixels.
 */
function calculateViewListMaxHeight(views: CategoryView[]): number {
  return views.reduce((acc, { label }) => {
    const listMaxHeight =
      (MAX_DISPLAYABLE_LIST_ITEMS.SINGLETON + 0.5) * VIEW_LIST_ITEM_HEIGHT +
      VIEW_LIST_SUBHEADER_HEIGHT * (label ? 1 : 0);
    return Math.max(acc, listMaxHeight);
  }, 0);
}

export default function FilterCategoryViews({
  categoryViews,
  onFilter,
}: Props): JSX.Element {
  // TODO(cc) filter category views to displayable views based on category view isDisabled, or OntologyCategoryTreeView has children.
  const viewListMaxHeight = calculateViewListMaxHeight(categoryViews);
  return (
    <ViewsMenu>
      {categoryViews.map((categoryView, i) => (
        <Fragment key={`${categoryView.key}${i}`}>
          {i !== 0 && <ViewDivider orientation="vertical" />}
          <FilterCategoryView
            categoryKey={categoryView.key}
            categoryView={categoryView}
            onFilter={onFilter}
            viewListMaxHeight={viewListMaxHeight}
          />
        </Fragment>
      ))}
    </ViewsMenu>
  );
}
