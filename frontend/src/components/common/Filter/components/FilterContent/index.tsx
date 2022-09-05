import { useEffect, useRef } from "react";
import {
  isOntologyCategoryView,
  isOntologyMultiPanelCategoryView,
  isSelectCategoryView,
} from "src/components/common/Filter";
import {
  CategoryView,
  OnFilterFn,
  OntologyMultiPanelCategoryView,
  SelectCategoryValueView,
} from "src/components/common/Filter/common/entities";
import FilterMenu from "src/components/common/Filter/components/FilterContent/components/FilterMenu";
import { MAX_DISPLAYABLE_MENU_ITEMS } from "src/components/common/Filter/components/FilterContent/components/FilterMenu/style";
import FilterRange from "src/components/common/Filter/components/FilterContent/components/FilterRange";
import FilterViews from "src/components/common/Filter/components/FilterContent/components/FilterViews";
import FilterMultiPanelCategoryView from "src/components/common/Filter/components/FilterContent/components/FilterViews/components/FilterMultiPanelCategoryView";
import {
  FilterSearchState,
  useFilterSearch,
} from "src/components/common/Filter/components/FilterSearch/common/useFilterSearch";
import { FilterContent as Content } from "./style";

interface Props {
  categoryView: CategoryView;
  onFilter: OnFilterFn;
}

export default function FilterContent({
  categoryView,
  onFilter,
}: Props): JSX.Element {
  const filterSearchState = useFilterSearch();
  const clientHeightRef = useRef<number>(0);
  const filterRef = useRef<HTMLDivElement>(null);
  const minHeight = clientHeightRef.current;

  useEffect(() => {
    if (filterRef.current) {
      clientHeightRef.current = filterRef.current.clientHeight;
    }
  }, []);

  return (
    <Content minHeight={minHeight} ref={filterRef}>
      {buildBasicFilterContent(categoryView, onFilter, filterSearchState)}
    </Content>
  );
}

/**
 * Build content model of basic filter depending on category type.
 * @param categoryView - View model of category to display.
 * @param onFilter - Function to execute on select of category value or remove of selected category value.
 * @param filterSearchState - Filter search value and corresponding functions to clear or update search value.
 * @returns Element representing content to display inside basic filter menu.
 */
function buildBasicFilterContent(
  categoryView: CategoryView,
  onFilter: OnFilterFn,
  filterSearchState: FilterSearchState
): JSX.Element {
  const { key } = categoryView;
  const { clearSearchValueFn, searchValue, setSearchValue } = filterSearchState;

  // Handle ontology categories.
  if (isOntologyCategoryView(categoryView)) {
    return (
      <FilterViews
        categoryFilterId={key}
        isSearchable={categoryView.isSearchable}
        isZerosVisible={categoryView.isZerosVisible}
        onFilter={onFilter}
        views={categoryView.views}
      />
    );
  }

  // Handle select categories
  if (isSelectCategoryView(categoryView)) {
    const { pinnedValues, unpinnedValues, values } = categoryView;
    return (
      <FilterMenu
        categoryFilterId={key}
        isMultiselect // Can possibly be single select with future filter types
        isSearchable={values.length > MAX_DISPLAYABLE_MENU_ITEMS}
        onFilter={onFilter}
        pinnedValues={filterCategoryValuesWithCount(pinnedValues)}
        searchValue={searchValue}
        setSearchValue={setSearchValue}
        unpinnedValues={filterCategoryValuesWithCount(unpinnedValues)}
        values={filterCategoryValuesWithCount(values)}
      />
    );
  }

  // Handle ontology multi panel categories
  if (isOntologyMultiPanelCategoryView(categoryView)) {
    return (
      <FilterMultiPanelCategoryView
        categoryView={categoryView}
        clearSearchValueFn={clearSearchValueFn}
        isSearchable={isFilterMultiPanelSearchable(categoryView)}
        onFilter={onFilter}
        searchValue={searchValue}
        setSearchValue={setSearchValue}
      />
    );
  }

  // Otherwise, handle range categories
  return <FilterRange categoryView={categoryView} onFilter={onFilter} />;
}

/**
 * Returns filtered category values where category count is greater than zero.
 * @param categoryValues - Category value view models for a given category.
 * @returns category values with a count.
 */
function filterCategoryValuesWithCount(
  categoryValues: SelectCategoryValueView[]
): SelectCategoryValueView[] {
  return categoryValues.filter(({ count }) => count > 0);
}

/**
 * Returns true if ontology multi panel is searchable.
 * @param categoryView - Ontology multi panel category view.
 * @returns true if ontology multi panel is searchable.
 */
function isFilterMultiPanelSearchable(
  categoryView: OntologyMultiPanelCategoryView
): boolean {
  let isSearchable = false;
  for (const panel of categoryView.panels) {
    const { views } = panel;
    if (views.length > MAX_DISPLAYABLE_MENU_ITEMS) {
      isSearchable = true;
      break;
    }
  }
  return isSearchable;
}
