import { Tooltip } from "@czi-sds/components";

import questionMarkIcon from "src/common/images/question-mark-icon.svg";
import { TooltipButton, StyledTooltip, StyledIconImage } from "./style";
import { Dispatch, ReactElement, SetStateAction } from "react";

interface Props {
  text: string | ReactElement;
  placement?:
    | "bottom-end"
    | "bottom-start"
    | "bottom"
    | "left-end"
    | "left-start"
    | "left"
    | "right-end"
    | "right-start"
    | "right"
    | "top-end"
    | "top-start"
    | "top";
  dark?: boolean;
  buttonDataTestId?: string;
  setTooltipContent: Dispatch<
    SetStateAction<{
      title: string;
      element: JSX.Element;
    } | null>
  >;
  title: string;
}
const HelpTooltip = ({
  text,
  placement = "right",
  dark,
  buttonDataTestId = "",
  setTooltipContent,
  title,
}: Props) => {
  return (
    <Tooltip
      sdsStyle={dark ? "dark" : "light"}
      placement={placement}
      width="wide"
      arrow
      title={<StyledTooltip>{text}</StyledTooltip>}
      slotProps={{
        tooltip: {
          style: dark ? { maxWidth: 550 } : {}, // Fixes SDS bug where "wide" property doesn't affect dark sdsStyle
        },
      }}
    >
      <TooltipButton
        data-testid={buttonDataTestId}
        sdsStyle="minimal"
        sdsType="secondary"
        isAllCaps={false}
      >
        <StyledIconImage
          onTouchEnd={() => {
            setTooltipContent({
              title: title,
              element: <StyledTooltip>{text}</StyledTooltip>,
            });
          }}
          src={questionMarkIcon}
        />
      </TooltipButton>
    </Tooltip>
  );
};

export default HelpTooltip;
