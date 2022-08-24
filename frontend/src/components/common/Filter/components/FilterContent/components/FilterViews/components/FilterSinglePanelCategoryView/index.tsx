import React from "react";
import { CategoryFilterId } from "src/common/hooks/useCategoryFilter";
import {
  OnFilterFn,
  OntologyPanelCategoryView,
  SelectCategoryValueView,
} from "src/components/common/Filter/common/entities";
import FilterCategoryViewPanel from "src/components/common/Filter/components/FilterContent/components/FilterViews/components/FilterCategoryViewPanel";
import FilterViewList from "src/components/common/Filter/components/FilterContent/components/FilterViews/components/FilterViewList";
import { ViewHeader } from "../FilterView/style";

interface Props {
  categoryKey: CategoryFilterId;
  categoryView: OntologyPanelCategoryView;
  onFilter: OnFilterFn;
  searchValue: string;
  viewListMaxHeight: number;
}

/**
 * Returns filtered select category views where category label includes search value.
 * @param views - Select category value views model for a given category.
 * @param searchValue - Search string to filter select category views.
 * @returns array of category views filtered by the given search value.
 */
function filterViews(
  views: SelectCategoryValueView[],
  searchValue: string
): SelectCategoryValueView[] {
  if (searchValue) {
    return views.filter(({ label }) =>
      label.toLowerCase().includes(searchValue)
    );
  }
  return views;
}

export default function FilterSinglePanelCategoryView({
  categoryKey,
  categoryView,
  onFilter,
  searchValue,
  viewListMaxHeight,
}: Props): JSX.Element {
  const { label, views } = categoryView;
  const filteredViews = filterViews(views, searchValue);
  return (
    <FilterCategoryViewPanel viewListMaxHeight={viewListMaxHeight}>
      <FilterViewList
        categoryKey={categoryKey}
        isZerosVisible={false}
        onFilter={onFilter}
        values={filteredViews}
        ViewHeader={label ? <ViewHeader>{label}</ViewHeader> : undefined}
      />
    </FilterCategoryViewPanel>
  );
}
