import React, { Fragment } from "react";
import {
  CategoryView,
  OnFilterFn,
} from "src/components/common/Filter/common/entities";
import { ViewDivider } from "src/components/common/Filter/components/FilterViews/components/FilterView/style";
import { ViewsMenu } from "src/components/common/Filter/components/FilterViews/style";
import FilterCategoryView from "../FilterCategoryView";

interface Props {
  categoryViews: CategoryView[];
  onFilter: OnFilterFn;
}

export default function FilterCategoryViews({
  categoryViews,
  onFilter,
}: Props): JSX.Element {
  // TODO(cc) filter category views to displayable views based on category view isDisabled, or OntologyCategoryTreeView has children.
  // TODO(cc) calculate the maximum possible pixel height of any view list.
  const viewListMaxHeight = 415; // TODO(cc)
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
