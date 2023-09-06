import { Tooltip } from "@czi-sds/components";

import questionMarkIcon from "src/common/images/question-mark-icon.svg";
import { TooltipButton, StyledTooltip, StyledIconImage } from "./style";
import { Dispatch, ReactElement, SetStateAction } from "react";

const getSlotProps = (dark?: boolean) => {
  return {
    tooltip: {
      style: dark ? { maxWidth: 550 } : {}, // Fixes SDS bug where "wide" property doesn't affect dark sdsStyle
    },
  };
};

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
  skinnyMode: boolean;
}
const HelpTooltip = ({
  text,
  placement = "right",
  dark,
  buttonDataTestId = "",
  setTooltipContent,
  title,
  skinnyMode = false,
}: Props) => {
  return (
    <Tooltip
      sdsStyle={dark ? "dark" : "light"}
      placement={placement}
      width="wide"
      arrow
      title={!skinnyMode && <StyledTooltip>{text}</StyledTooltip>}
      slotProps={getSlotProps(dark)}
    >
      <TooltipButton
        data-testid={buttonDataTestId}
        sdsStyle="minimal"
        sdsType="secondary"
        isAllCaps={false}
      >
        <StyledIconImage
          onClick={() => {
            if (skinnyMode) {
              setTooltipContent({
                title: title,
                element: <StyledTooltip>{text}</StyledTooltip>,
              });
            }
          }}
          src={questionMarkIcon}
        />
      </TooltipButton>
    </Tooltip>
  );
};

export default HelpTooltip;
