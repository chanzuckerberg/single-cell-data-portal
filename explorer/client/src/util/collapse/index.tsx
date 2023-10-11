/* Core dependencies */
import { Colors, Icon } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import React, { FC, ReactNodeArray, useState } from "react";

/* Styles */
import * as globals from "../../globals";

interface Props {
  children: ReactNodeArray;
  isOpen?: boolean;
}

/*
 Expands and collapses content, executed by onClick or onKeyPress of the collapse target.
 */
const Collapse: FC<Props> = ({ children, isOpen = true }): JSX.Element => {
  const [label, content] = children;
  const [isExpanded, setIsExpanded] = useState(isOpen);
  return (
    <>
      <span
        onClick={() => setIsExpanded((expanded) => !expanded)}
        onKeyPress={() => setIsExpanded((expanded) => !expanded)}
        role="menuitem"
        style={{
          alignItems: "center",
          color: Colors.BLACK,
          cursor: "pointer",
          display: "flex",
          flexBasis: "100%",
          lineHeight: "16px",
          margin: "8px 0",
          width: "fit-content",
        }}
        tabIndex={0}
      >
        <span
          style={{ fontWeight: globals.bold, marginRight: 4, padding: "3px 0" }}
        >
          {label}
        </span>
        {isExpanded ? (
          <Icon icon={IconNames.CHEVRON_DOWN} />
        ) : (
          <Icon icon={IconNames.CHEVRON_RIGHT} />
        )}
      </span>
      {isExpanded ? content : null}
    </>
  );
};

export default Collapse;
