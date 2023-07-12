import { ListSubheader } from "@czi-sds/components";
import React, { useEffect, useRef, useState } from "react";
import { useResizeObserver } from "src/common/hooks/useResizeObserver";
import {
  CATEGORY_FILTER_ID,
  OnFilterFn,
  OntologyCategoryTreeNodeView,
} from "src/components/common/Filter/common/entities";
import {
  ViewDivider,
  ViewPanel,
  ViewPanelScroll,
} from "src/components/common/Filter/components/FilterContent/components/FilterViews/components/FilterView/style";
import FilterViewList from "src/components/common/Filter/components/FilterContent/components/FilterViews/components/FilterViewList";
import FilterSearch from "src/components/common/Filter/components/FilterSearch";
import { useFilterSearch } from "src/components/common/Filter/components/FilterSearch/common/useFilterSearch";
import { ADDITIONAL_MENU_WIDTH } from "src/components/common/Filter/components/FilterContent/components/FilterMenu";

export const MAX_DISPLAYABLE_LIST_ITEMS = {
  NON_SINGLETON: 15,
  SINGLETON: 9,
};

interface Props {
  categoryFilterId: CATEGORY_FILTER_ID;
  label?: string;
  isSearchable: boolean;
  isZerosVisible: boolean;
  onFilter: OnFilterFn;
  showViewDivider: boolean;
  values: OntologyCategoryTreeNodeView[];
  viewListMaxHeight: number;
}

export default function FilterView({
  categoryFilterId,
  label,
  isSearchable,
  isZerosVisible,
  onFilter,
  showViewDivider,
  values,
  viewListMaxHeight,
}: Props): JSX.Element {
  const { searchValue, setSearchValue } = useFilterSearch();
  const panelRef = useRef<HTMLDivElement>(null);
  const listContainerRef = useRef<HTMLDivElement>(null);
  const listContainerRect = useResizeObserver(listContainerRef);
  const [panelScrollable, setPanelScrollable] = useState(false);
  const [panelWidth, setPanelWidth] = useState<number>(0);
  const { scrollHeight: listScrollHeight } = listContainerRect || {};
  const filteredValues = filterViewValues(values, searchValue);

  // Calculate and set a min width on view panel to prevent width resizing
  // (derived from a change in list item selected state font weight).
  useEffect(() => {
    if (panelRef.current) {
      setPanelWidth(panelRef.current?.clientWidth + ADDITIONAL_MENU_WIDTH);
    }
  }, []);

  // Set scrollable state on change of filter view list scroll height.
  useEffect(() => {
    if (listScrollHeight) {
      setPanelScrollable(listScrollHeight > viewListMaxHeight);
    }
  }, [listScrollHeight, viewListMaxHeight]);

  return (
    <>
      {showViewDivider && <ViewDivider orientation="vertical" />}
      <ViewPanel panelWidth={panelWidth} ref={panelRef}>
        {/* Optional search bar */}
        {isSearchable && (
          <FilterSearch
            searchValue={searchValue}
            setSearchValue={setSearchValue}
          />
        )}
        <ViewPanelScroll
          maxHeight={viewListMaxHeight}
          ref={listContainerRef}
          scrollable={panelScrollable}
        >
          <FilterViewList
            categoryFilterId={categoryFilterId}
            isZerosVisible={isZerosVisible}
            onFilter={onFilter}
            values={filteredValues}
            ViewHeader={
              label ? <ListSubheader>{label}</ListSubheader> : undefined
            }
          />
        </ViewPanelScroll>
      </ViewPanel>
    </>
  );
}

/**
 * Returns filtered ontology tree view values where value label includes search value.
 * @param values - Ontology tree view values.
 * @param searchValue - Search string that filters category values.
 * @returns array of ontology tree view values filtered by the given search value.
 */
function filterViewValues(
  values: OntologyCategoryTreeNodeView[] = [],
  searchValue: string
): OntologyCategoryTreeNodeView[] {
  if (searchValue) {
    return [...values].reduce((acc, value) => {
      // Clone the node and any corresponding children.
      const nodeValue = { ...value };
      const nodeChildren = nodeValue.children;
      if (nodeChildren) {
        // Update the node with filtered children by the search value (recursive).
        nodeValue.children = filterViewValues(nodeChildren, searchValue);
        if (nodeValue.children.length > 0) {
          // Accumulate the node; filtered children include the search value.
          acc.push(nodeValue);
          return acc;
        }
      }
      if (isSearchHit(nodeValue.label, searchValue)) {
        // Otherwise, if the node label includes the search value, delete the children property - filtered children
        // have returned an empty array - and accumulate the modified node.
        delete nodeValue.children;
        acc.push(nodeValue);
        return acc;
      }
      return acc;
    }, [] as OntologyCategoryTreeNodeView[]);
  }
  return values;
}

/**
 * Returns true if the value includes the search string.
 * @param value - View list value.
 * @param searchValue - View search value.
 * @returns true if the value includes the search string.
 */
function isSearchHit(value: string, searchValue: string): boolean {
  return value.trim().toLowerCase().includes(searchValue.toLowerCase());
}
