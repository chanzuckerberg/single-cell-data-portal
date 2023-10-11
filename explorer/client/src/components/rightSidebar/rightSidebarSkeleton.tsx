/* Core dependencies */
import { SKELETON } from "@blueprintjs/core/lib/esnext/common/classes";
import React from "react";

/* Styles */
import { STYLE_RIGHT_SIDEBAR } from ".";

/**
 * Skeleton of right side bar, to be displayed during data load.
 * @returns Markup displaying right side bar skeleton.
 */
function RightSidebarSkeleton(): JSX.Element {
  return (
    <div style={{ ...STYLE_RIGHT_SIDEBAR, top: -2 }}>
      {/* Quick gene search */}
      {/* Gene menu */}
      <div
        style={{
          height: 21,
          marginBottom: 10,
          position: "relative",
          width: 65,
        }}
        className={SKELETON}
      />
      {/* Gene search <input/> */}
      <div style={{ height: 30, marginBottom: 16 }} className={SKELETON} />
      {/* Create new gene set button */}
      <div style={{ display: "flex", justifyContent: "space-between" }}>
        {/* Gene sets menu */}
        <div
          style={{
            height: 30,
            width: 93,
          }}
          className={SKELETON}
        />
        {/* Create new button */}
        <div
          style={{
            height: 30,
            width: 82,
          }}
          className={SKELETON}
        />
      </div>
    </div>
  );
}

export default RightSidebarSkeleton;
