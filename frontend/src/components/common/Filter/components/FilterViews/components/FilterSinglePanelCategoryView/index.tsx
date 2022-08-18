import React from "react";
import { CategoryFilterId } from "src/common/hooks/useCategoryFilter";
import {
  OnFilterFn,
  OntologyPanelCategoryView,
} from "src/components/common/Filter/common/entities";
import FilterCategoryViewPanel from "src/components/common/Filter/components/FilterViews/components/FilterCategoryViewPanel";
import FilterViewList from "src/components/common/Filter/components/FilterViews/components/FilterViewList";
import { ViewHeader } from "../FilterView/style";

interface Props {
  categoryKey: CategoryFilterId;
  categoryView: OntologyPanelCategoryView;
  onFilter: OnFilterFn;
  viewListMaxHeight: number;
}

export default function FilterSinglePanelCategoryView({
  categoryKey,
  categoryView,
  onFilter,
  viewListMaxHeight,
}: Props): JSX.Element {
  const { label, views } = categoryView;
  return (
    <FilterCategoryViewPanel viewListMaxHeight={viewListMaxHeight}>
      <FilterViewList
        categoryKey={categoryKey}
        isZerosVisible={false}
        onFilter={onFilter}
        values={views}
        ViewHeader={
          label ? <ViewHeader disableSticky>{label}</ViewHeader> : undefined
        }
      />
    </FilterCategoryViewPanel>
  );
}
