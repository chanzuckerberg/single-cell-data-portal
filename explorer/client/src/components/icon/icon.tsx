/* core dependencies */
import React from "react";
import { Classes, IconSize } from "@blueprintjs/core";

/* app dependencies */
import { IconName } from "./iconName";
import { IconSvgPaths } from "./iconSvgPaths";

interface Props {
  icon: IconName;
}

/**
 * Renders custom icons outside of standard Blueprint icon set. Implementation based on Blueprint's Icon component.
 */
function Icon(props: Props): JSX.Element | null {
  const { icon } = props;

  /*
   Determine the SVG path configuration for the given icon name.
   @param iconName - Name of the icon to render.
   @returns Configured path element for the given icon name.
   */
  const renderSvgPaths = (iconName: IconName): JSX.Element[] | null => {
    const paths = IconSvgPaths[iconName];
    if (!paths) {
      return null;
    }
    return paths.map(({ d, fill }) => (
      <path key={d} d={d} fill={fill} fillRule="evenodd" />
    ));
  };

  const iconSize = IconSize.STANDARD;
  const viewBox = `0 0 ${iconSize} ${iconSize}`;
  const paths = renderSvgPaths(icon);

  return paths ? (
    <span className={Classes.ICON}>
      <svg
        fill="none"
        data-icon={icon}
        height={iconSize}
        viewBox={viewBox}
        width={iconSize}
      >
        {paths}
      </svg>
    </span>
  ) : null;
}

export default Icon;
