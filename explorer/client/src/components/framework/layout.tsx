import React, { Children, useState } from "react";
import * as globals from "../../globals";

interface Props {
  datasetMetadataError: string | null;
  renderGraph?: (viewport: HTMLDivElement) => React.ReactNode;
}

/*
  Layout - this react component contains all the layout style and logic for the application once it has loaded.

  The layout is based on CSS grid: the left and right sidebars have fixed widths, the graph in the middle takes the
  remaining space.

  Note, the renderGraph child is a function rather than a fully-instantiated element because the middle pane of the
  app is dynamically-sized. It must have access to the containing viewport in order to know how large the graph
  should be.
*/

const Layout: React.FC<Props> = (props) => {
  const [viewportRef, setViewportRef] = useState<HTMLDivElement | null>(null);

  const { children, datasetMetadataError, renderGraph } = props;
  const [leftSidebar, rightSidebar] = Children.toArray(children);
  let graphComponent = null;
  if (viewportRef && renderGraph) {
    graphComponent = renderGraph(viewportRef);
  }
  return (
    <div
      style={{
        display: "grid",
        paddingTop: !datasetMetadataError ? globals.HEADER_HEIGHT_PX : 0,
        gridTemplateColumns: `
        [left-sidebar-start] ${globals.leftSidebarWidth + 1}px
        [left-sidebar-end graph-start] auto
        [graph-end right-sidebar-start] ${
          globals.rightSidebarWidth + 1
        }px [right-sidebar-end]
      `,
        gridTemplateRows: "[top] auto [bottom]",
        gridTemplateAreas: "left-sidebar | graph | right-sidebar",
        columnGap: "0px",
        justifyItems: "stretch",
        alignItems: "stretch",
        height: "inherit",
        width: "inherit",
        position: "relative",
        top: 0,
        left: 0,
        minWidth: "1240px",
      }}
    >
      <div
        style={{
          gridArea: "top / left-sidebar-start / bottom / left-sidebar-end",
          position: "relative",
          height: "inherit",
          overflowY: "auto",
        }}
      >
        {leftSidebar}
      </div>
      <div
        style={{
          zIndex: 0,
          gridArea: "top / graph-start / bottom / graph-end",
          position: "relative",
          height: "inherit",
        }}
        ref={(ref) => {
          setViewportRef(ref);
        }}
      >
        {graphComponent}
      </div>
      <div
        style={{
          gridArea: "top / right-sidebar-start / bottom / right-sidebar-end",
          position: "relative",
          height: "inherit",
          overflowY: "auto",
        }}
      >
        {/* The below conditional is required because the right sidebar initializes as function for some reason...*/}
        {!(rightSidebar instanceof Function) && rightSidebar}
      </div>
    </div>
  );
};

Layout.defaultProps = {
  renderGraph: undefined,
};
export default Layout;
