// TODO(cc) deprecated FilterCategoryView.
import React from "react";
import { CategoryFilterId } from "src/common/hooks/useCategoryFilter";
import {
  isOntologyCategoryView,
  isRangeCategoryView,
  isSelectCategoryView,
} from "src/components/common/Filter";
import {
  CategoryView,
  OnFilterFn,
} from "src/components/common/Filter/common/entities";
import FilterRange from "src/components/common/Filter/components/FilterRange";
import FilterCategoryViewPanel from "src/components/common/Filter/components/FilterViews/components/FilterCategoryViewPanel";
import FilterViewList from "src/components/common/Filter/components/FilterViews/components/FilterViewList";
import { ViewHeader } from "../FilterView/style";

interface Props {
  categoryKey: CategoryFilterId;
  categoryView: CategoryView;
  onFilter: OnFilterFn;
  viewListMaxHeight: number;
}

export default function FilterCategoryView({
  categoryKey,
  categoryView,
  onFilter,
  viewListMaxHeight,
}: Props): JSX.Element {
  const { label } = categoryView;
  return (
    <FilterCategoryViewPanel viewListMaxHeight={viewListMaxHeight}>
      {/* Select category view */}
      {isSelectCategoryView(categoryView) && (
        <FilterViewList
          categoryKey={categoryKey}
          isZerosVisible={false}
          onFilter={onFilter}
          values={categoryView.values}
          ViewHeader={
            label ? <ViewHeader disableSticky>{label}</ViewHeader> : undefined
          }
        />
      )}
      {/* Ontology category view */}
      {isOntologyCategoryView(categoryView) && (
        <FilterViewList
          categoryKey={categoryKey}
          isZerosVisible={categoryView.isZerosVisible}
          onFilter={onFilter}
          values={categoryView.views[0].children}
          ViewHeader={
            label ? <ViewHeader disableSticky>{label}</ViewHeader> : undefined
          }
        />
      )}
      {/* Range category view */}
      {isRangeCategoryView(categoryView) && (
        <FilterRange categoryView={categoryView} onFilter={onFilter} />
      )}
    </FilterCategoryViewPanel>
  );
}
