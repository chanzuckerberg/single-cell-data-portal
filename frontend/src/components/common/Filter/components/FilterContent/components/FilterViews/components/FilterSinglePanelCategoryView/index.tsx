import { ListSubheader } from "@czi-sds/components";
import React from "react";
import {
  CATEGORY_FILTER_ID,
  OnFilterFn,
  OntologyPanelCategoryView,
  SelectCategoryValueView,
} from "src/components/common/Filter/common/entities";
import FilterCategoryViewPanel from "src/components/common/Filter/components/FilterContent/components/FilterViews/components/FilterCategoryViewPanel";
import FilterViewList from "src/components/common/Filter/components/FilterContent/components/FilterViews/components/FilterViewList";

interface Props {
  categoryFilterId: CATEGORY_FILTER_ID;
  categoryView: OntologyPanelCategoryView;
  onFilter: OnFilterFn;
  searchValue: string;
  viewListMaxHeight: number;
}

/**
 * Returns filtered select category views where category label includes search value.
 * If there is no search value, only "visible" views are returned.
 * @param views - Select category value views model for a given category.
 * @param searchValue - Search string to filter select category views.
 * @returns array of category views filtered by the given search value.
 */
function filterViews(
  views: SelectCategoryValueView[],
  searchValue: string
): SelectCategoryValueView[] {
  if (searchValue) {
    const lowerCaseSearchValue = searchValue.toLowerCase();
    return views.filter(({ label }) =>
      label.toLowerCase().includes(lowerCaseSearchValue)
    );
  } else {
    return views.filter((view) => view.visible);
  }
}

export default function FilterSinglePanelCategoryView({
  categoryFilterId,
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
        categoryFilterId={categoryFilterId}
        isZerosVisible={false}
        onFilter={onFilter}
        values={filteredViews}
        ViewHeader={label ? <ListSubheader>{label}</ListSubheader> : undefined}
      />
    </FilterCategoryViewPanel>
  );
}
