/* Core dependencies */
import { SKELETON } from "@blueprintjs/core/lib/esnext/common/classes";
import React from "react";

/* App dependencies */
import Controls from "../controls";
import Layout from "./layout";
import LeftSidebarSkeleton from "../leftSidebar/leftSidebarSkeleton";
import RightSidebarSkeleton from "../rightSidebar/rightSidebarSkeleton";

/* Padding between dataset selector and menubar */
const PADDING_CONTROLS = 10;

/* Width of menubar */
const WIDTH_MENUBAR = 482;

/**
 * Skeleton layout component displayed when in loading state.
 * @returns Markup displaying skeleton.
 */
function LayoutSkeleton(): JSX.Element {
  return (
    <Layout datasetMetadataError={null}>
      <LeftSidebarSkeleton />
      {() => (
        <Controls>
          <div
            style={{
              height: 30,
              position: "relative",
              top: 8,
              width: `calc(100% - ${WIDTH_MENUBAR}px - ${PADDING_CONTROLS}px)`,
            }}
            className={SKELETON}
          />
          <div
            style={{
              height: 30,
              position: "relative",
              top: 8,
              width: WIDTH_MENUBAR,
            }}
            className={SKELETON}
          />
        </Controls>
      )}
      <RightSidebarSkeleton />
    </Layout>
  );
}

export default LayoutSkeleton;
