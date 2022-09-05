import React, { Fragment } from "react";
import {
  CategoryValueId,
  CATEGORY_FILTER_ID,
  MultiPanelOntologyCategoryView,
  OnFilterFn,
  OntologyPanelCategoryView,
  Range,
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
import {
  ClearSearchValueFn,
  SetSearchValueFn,
} from "src/components/common/Filter/components/FilterSearch/common/useFilterSearch";

interface Props {
  categoryView: MultiPanelOntologyCategoryView;
  clearSearchValueFn: ClearSearchValueFn;
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
  clearSearchValueFn,
  isSearchable,
  onFilter,
  searchValue,
  setSearchValue,
}: Props): JSX.Element {
  const { categoryFilterId, panels } = categoryView;
  const viewListMaxHeight = calculateViewListMaxHeight(panels);
  // Update onFilter function with clear search value function.
  const onFilterWithClearSearch = (
    categoryFilterId: CATEGORY_FILTER_ID,
    key: CategoryValueId | null, // null for ranges.
    value: CategoryValueId | Range
  ) => {
    onFilter(categoryFilterId, key, value);
    clearSearchValueFn();
  };
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
        {panels.map((ontologyCategoryView, i) => {
          const { isSearchMultiselect } = ontologyCategoryView;
          const onFilterFn = isSearchMultiselect
            ? onFilter
            : onFilterWithClearSearch;
          return (
            <Fragment key={`${ontologyCategoryView.label}${i}`}>
              {i !== 0 && <ViewDivider orientation="vertical" />}
              <FilterSinglePanelCategoryView
                categoryFilterId={categoryFilterId}
                categoryView={ontologyCategoryView}
                onFilter={onFilterFn}
                searchValue={searchValue}
                viewListMaxHeight={viewListMaxHeight}
              />
            </Fragment>
          );
        })}
      </ViewsMenu>
    </Views>
  );
}
