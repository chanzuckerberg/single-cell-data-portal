import { Tooltip } from "@czi-sds/components";
import { StyledTooltip, TooltipContent, TooltipTrigger } from "../../../style";
import { ReactNode } from "react";

export function FilterFooterTooltip({
  children,
  content,
}: {
  children: ReactNode;
  content: ReactNode;
}) {
  return (
    <Tooltip
      sdsStyle="dark"
      placement="top"
      arrow
      slotProps={{
        tooltip: {
          style: {
            maxWidth: 395, // Override the max-width specification for dark sdsStyle.
          },
        },
      }}
      title={
        <StyledTooltip>
          <TooltipContent>{content}</TooltipContent>
        </StyledTooltip>
      }
    >
      <TooltipTrigger>{children}</TooltipTrigger>
    </Tooltip>
  );
}
