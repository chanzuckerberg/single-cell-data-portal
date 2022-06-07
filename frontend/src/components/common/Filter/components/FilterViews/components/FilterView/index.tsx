import React, { useEffect, useRef, useState } from "react";
import { useResizeObserver } from "src/common/hooks/useResizeObserver";
import {
  CATEGORY_KEY,
  OnFilterFn,
  OntologyCategoryTreeNodeView,
  OnUpdateSearchValueFn,
} from "src/components/common/Filter/common/entities";
import {
  ViewDivider,
  ViewHeader,
  ViewPanel,
  ViewPanelScroll,
} from "src/components/common/Filter/components/FilterViews/components/FilterView/style";
import FilterViewList from "src/components/common/Filter/components/FilterViews/components/FilterViewList";
import FilterViewSearch from "src/components/common/Filter/components/FilterViews/components/FilterViewSearch";

const ADDITIONAL_PANEL_WIDTH = 8;

export const MAX_DISPLAYABLE_LIST_ITEMS = {
  NON_SINGLETON: 15,
  SINGLETON: 9,
};

interface Props {
  categoryKey: CATEGORY_KEY;
  label?: string;
  isSearchable: boolean;
  isZerosVisible: boolean;
  onFilter: OnFilterFn;
  onUpdateSearchValue: OnUpdateSearchValueFn;
  showViewDivider: boolean;
  values: OntologyCategoryTreeNodeView[];
  viewListMaxHeight: number;
}

export default function FilterView({
  categoryKey,
  label,
  isSearchable,
  isZerosVisible,
  onFilter,
  onUpdateSearchValue,
  showViewDivider,
  values,
  viewListMaxHeight,
}: Props): JSX.Element {
  const listContainerRef = useRef<HTMLDivElement>(null);
  const listContainerRect = useResizeObserver(listContainerRef);
  const [panelScrollable, setPanelScrollable] = useState(false);
  const [panelWidth, setPanelWidth] = useState<number>(0);
  const [searchValue, setSearchValue] = useState<string>("");
  const { scrollHeight: listScrollHeight } = listContainerRect || {};
  const filteredValues = filterViewValues(values, searchValue);

  // Calculate and set a min width on view list to prevent width resizing
  // (derived from a change in list item selected state font weight).
  useEffect(() => {
    if (listContainerRef.current) {
      setPanelWidth(
        listContainerRef.current?.clientWidth + ADDITIONAL_PANEL_WIDTH
      );
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
      <ViewPanel panelWidth={panelWidth}>
        {/* Optional search bar */}
        {isSearchable && (
          <FilterViewSearch
            onUpdateSearchValue={onUpdateSearchValue}
            setSearchValue={setSearchValue}
          />
        )}
        <ViewPanelScroll
          maxHeight={viewListMaxHeight}
          ref={listContainerRef}
          scrollable={panelScrollable}
        >
          <FilterViewList
            categoryKey={categoryKey}
            isZerosVisible={isZerosVisible}
            onFilter={onFilter}
            values={filteredValues}
            ViewHeader={
              label ? <ViewHeader disableSticky>{label}</ViewHeader> : undefined
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
