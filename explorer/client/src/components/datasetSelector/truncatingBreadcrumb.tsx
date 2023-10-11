/* Core dependencies */
import { Breadcrumb } from "@blueprintjs/core";
import React from "react";

/* App dependencies */
import { TruncatingBreadcrumbProps } from "./truncatingBreadcrumbs";

// Styles
// @ts-expect-error --- TODO fix import
import styles from "./datasetSelector.css";

interface Props {
  children?: React.ReactNode;
  item: TruncatingBreadcrumbProps;
}

/**
 * Breadcrumb displayed in truncating breadcrumbs.
 */
const TruncatingBreadcrumb = React.memo<Props>(
  ({ item }): JSX.Element => (
    <Breadcrumb
      href={item.href}
      onClick={item.disabled ? undefined : item.onClick}
      className={
        item.disabled
          ? styles.datasetSelectorBreadcrumbDisabled
          : styles.datasetSelectorBreadcrumb
      }
    >
      {item.text}
      {item.icon}
    </Breadcrumb>
  ),
);
export default TruncatingBreadcrumb;
