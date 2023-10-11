import { SKELETON } from "@blueprintjs/core/lib/esnext/common/classes";
import React from "react";

interface StillLoadingProps {
  height?: number;
}
/**
 * Render a loading indicator for the field.
 */
const StillLoading = ({ height = 211 }: StillLoadingProps): JSX.Element => (
  <div style={{ height }} className={SKELETON} />
);
StillLoading.defaultProps = {
  height: 211,
};

export default StillLoading;
