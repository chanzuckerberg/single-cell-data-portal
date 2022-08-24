import React, { ReactNode, useEffect, useRef, useState } from "react";
import { useResizeObserver } from "src/common/hooks/useResizeObserver";
import {
  CategoryViewPanel,
  ViewPanelScroll,
} from "src/components/common/Filter/components/FilterContent/components/FilterViews/components/FilterView/style";

interface Props {
  children: ReactNode;
  viewListMaxHeight: number;
}

export default function FilterCategoryViewPanel({
  children,
  viewListMaxHeight,
}: Props): JSX.Element {
  const listContainerRef = useRef<HTMLDivElement>(null);
  const listContainerRect = useResizeObserver(listContainerRef);
  const [panelScrollable, setPanelScrollable] = useState(false);
  const { scrollHeight: listScrollHeight } = listContainerRect || {};

  // Set scrollable state on change of filter view list scroll height.
  useEffect(() => {
    if (listScrollHeight) {
      setPanelScrollable(listScrollHeight > viewListMaxHeight);
    }
  }, [listScrollHeight, viewListMaxHeight]);

  return (
    <CategoryViewPanel>
      <ViewPanelScroll
        maxHeight={viewListMaxHeight}
        ref={listContainerRef}
        scrollable={panelScrollable}
      >
        {children}
      </ViewPanelScroll>
    </CategoryViewPanel>
  );
}
