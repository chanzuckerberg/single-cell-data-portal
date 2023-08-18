import { Tooltip } from "@czi-sds/components";
import { ReactElement, useMemo } from "react";
import { categoryTooltipCss } from "src/components/common/Filter/components/FilterLabel/components/FilterLabelTooltip/style";

interface Props {
  children: ReactElement;
  tooltip?: string;
}

export default function FilterLabelTooltip({
  children,
  tooltip,
}: Props): JSX.Element {
  const tooltipClasses = useMemo(() => ({ tooltip: categoryTooltipCss }), []);
  return tooltip ? (
    <Tooltip
      arrow
      classes={tooltipClasses}
      placement="top"
      sdsStyle="dark"
      title={tooltip}
    >
      {children}
    </Tooltip>
  ) : (
    children
  );
}
