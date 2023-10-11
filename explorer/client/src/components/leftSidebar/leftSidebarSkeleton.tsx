/* Core dependencies */
import { SKELETON } from "@blueprintjs/core/lib/esnext/common/classes";
import React, { CSSProperties } from "react";

/* Styles */
import { STYLE_LEFT_SIDEBAR } from ".";
import StillLoading from "../brushableHistogram/loading";
import { StillLoading as CategoryLoading } from "../categorical/category";

const STYLE_SUPER_CATEGORY: CSSProperties = {
  height: 22,
  margin: "8px 0",
  width: 160,
};

/**
 * Skeleton of left side bar, to be displayed during data load.
 * @returns Markup displaying left side bar skeleton.
 */
function LeftSidebarSkeleton(): JSX.Element {
  return (
    <div style={STYLE_LEFT_SIDEBAR}>
      {/* Categorical */}
      <div style={{ padding: 10 }}>
        <div style={STYLE_SUPER_CATEGORY} className={SKELETON} />
        {[...Array(10)].map((_, i) => (
          /* eslint-disable-next-line react/no-array-index-key -- array order won't change */
          <CategoryLoading key={i} />
        ))}
      </div>
      {/* Continuous */}
      <div>
        <div
          style={{ ...STYLE_SUPER_CATEGORY, marginLeft: 10 }}
          className={SKELETON}
        />
        {[...Array(2)].map((_, i) => (
          /* eslint-disable-next-line react/no-array-index-key -- array order won't change */
          <StillLoading key={i} />
        ))}
      </div>
    </div>
  );
}

export default LeftSidebarSkeleton;
