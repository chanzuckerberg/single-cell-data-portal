import React, { Fragment } from "react";
import {
  OnFilterFn,
  OntologyMultiPanelCategoryView,
  OntologyPanelCategoryView,
} from "src/components/common/Filter/common/entities";
import {
  Search as FilterSearch,
  Views,
} from "src/components/common/Filter/components/FilterContent/components/FilterViews/components/FilterMultiPanelCategoryView/style";
import FilterSinglePanelCategoryView from "src/components/common/Filter/components/FilterContent/components/FilterViews/components/FilterSinglePanelCategoryView";
import { MAX_DISPLAYABLE_LIST_ITEMS } from "src/components/common/Filter/components/FilterContent/components/FilterViews/components/FilterView";
import {
  ViewDivider,
  VIEW_LIST_ITEM_HEIGHT,
  VIEW_LIST_SUBHEADER_HEIGHT,
} from "src/components/common/Filter/components/FilterContent/components/FilterViews/components/FilterView/style";
import { ViewsMenu } from "src/components/common/Filter/components/FilterContent/components/FilterViews/style";
import { SetSearchValueFn } from "src/components/common/Filter/components/FilterSearch/common/useFilterSearch";

interface Props {
  categoryView: OntologyMultiPanelCategoryView;
  isSearchable: boolean;
  onFilter: OnFilterFn;
  searchValue: string;
  setSearchValue: SetSearchValueFn;
}

/**
 * Returns the maximum possible pixel height of any view list.
 * @param views - Array of views within an ontology-aware filter.
 * @return maximum possible view list height in pixels.
 */
function calculateViewListMaxHeight(
  views: OntologyPanelCategoryView[]
): number {
  return views.reduce((acc, { label }) => {
    const listMaxHeight =
      (MAX_DISPLAYABLE_LIST_ITEMS.SINGLETON + 0.5) * VIEW_LIST_ITEM_HEIGHT +
      VIEW_LIST_SUBHEADER_HEIGHT * (label ? 1 : 0);
    return Math.max(acc, listMaxHeight);
  }, 0);
}

export default function FilterMultiPanelCategoryView({
  categoryView,
  isSearchable,
  onFilter,
  searchValue,
  setSearchValue,
}: Props): JSX.Element {
  const { key, panels } = categoryView;
  const viewListMaxHeight = calculateViewListMaxHeight(panels);
  return (
    <Views>
      {/* Optional search bar */}
      {isSearchable && (
        <FilterSearch
          searchValue={searchValue}
          setSearchValue={setSearchValue}
        />
      )}
      <ViewsMenu>
        {panels.map((ontologyCategoryView, i) => (
          <Fragment key={`${ontologyCategoryView.label}${i}`}>
            {i !== 0 && <ViewDivider orientation="vertical" />}
            <FilterSinglePanelCategoryView
              categoryKey={key}
              categoryView={ontologyCategoryView}
              onFilter={onFilter}
              searchValue={searchValue}
              viewListMaxHeight={viewListMaxHeight}
            />
          </Fragment>
        ))}
      </ViewsMenu>
    </Views>
  );
}
