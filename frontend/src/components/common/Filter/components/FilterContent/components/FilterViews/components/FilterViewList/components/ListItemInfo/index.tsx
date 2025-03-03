import { Tooltip } from "@czi-sds/components";
import { StyledTooltip, TooltipTrigger } from "./style";
import { ReactNode } from "react";

export function ListItemInfo({
  children,
  content,
}: {
  children: ReactNode;
  content: ReactNode;
}) {
  return (
    <TooltipTrigger>
      <Tooltip
        sdsStyle="dark"
        placement="top"
        arrow
        slotProps={{
          tooltip: {
            style: {
              maxWidth: 245, // Override the max-width specification for dark sdsStyle.
            },
          },
        }}
        title={<StyledTooltip>{content}</StyledTooltip>}
      >
        <span>{children}</span>
      </Tooltip>
    </TooltipTrigger>
  );
}
