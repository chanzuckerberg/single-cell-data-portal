import React, { ReactNode, useEffect, useRef, useState } from "react";
import { useResizeObserver } from "src/common/hooks/useResizeObserver";
import {
  ViewPanel,
  ViewPanelScroll,
} from "src/components/common/Filter/components/FilterViews/components/FilterView/style";

const ADDITIONAL_PANEL_WIDTH = 8; // TODO(cc) review duplication with FilterView

interface Props {
  children: ReactNode;
  viewListMaxHeight: number;
}

export default function FilterCategoryViewPanel({
  children,
  viewListMaxHeight,
}: Props): JSX.Element {
  const panelRef = useRef<HTMLDivElement>(null);
  const listContainerRef = useRef<HTMLDivElement>(null);
  const listContainerRect = useResizeObserver(listContainerRef);
  const [panelScrollable, setPanelScrollable] = useState(false);
  const [panelWidth, setPanelWidth] = useState<number>(0);
  const { scrollHeight: listScrollHeight } = listContainerRect || {};

  // Calculate and set a min width on view panel to prevent width resizing
  // (derived from a change in list item selected state font weight).
  useEffect(() => {
    if (panelRef.current) {
      setPanelWidth(panelRef.current?.clientWidth + ADDITIONAL_PANEL_WIDTH);
    }
  }, []);

  // Set scrollable state on change of filter view list scroll height.
  useEffect(() => {
    if (listScrollHeight) {
      setPanelScrollable(listScrollHeight > viewListMaxHeight);
    }
  }, [listScrollHeight, viewListMaxHeight]);

  return (
    <ViewPanel panelWidth={panelWidth} ref={panelRef}>
      <ViewPanelScroll
        maxHeight={viewListMaxHeight}
        ref={listContainerRef}
        scrollable={panelScrollable}
      >
        {children}
      </ViewPanelScroll>
    </ViewPanel>
  );
}
