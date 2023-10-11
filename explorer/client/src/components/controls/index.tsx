/* Core dependencies */
import React from "react";

interface Props {
  children: React.ReactNode;
  bottom?: number;
}

/**
 * Controls component for positioning graph controls.
 * @returns Markup displaying children positioned within the graph grid template area.
 */
function Controls(props: Props): JSX.Element {
  const { children, bottom } = props;
  return (
    <div
      style={{
        display: "flex",
        justifyContent: "space-between",
        left: 8,
        position: "absolute",
        right: 8,
        zIndex: 3,
        ...(bottom !== undefined ? { bottom } : { top: 0 }),
      }}
    >
      {children}
    </div>
  );
}

Controls.defaultProps = {
  bottom: undefined,
};

export default Controls;
