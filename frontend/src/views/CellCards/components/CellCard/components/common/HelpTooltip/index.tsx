import { Tooltip } from "@czi-sds/components";

import questionMarkIcon from "src/common/images/question-mark-icon.svg";
import { TooltipButton, StyledTooltip, StyledIconImage } from "./style";
import { ReactElement } from "react";

interface Props {
  text: string | ReactElement;
  placement?: "top" | "bottom" | "left" | "right";
  dark?: boolean;
  buttonDataTestId?: string;
}
const HelpTooltip = ({
  text,
  placement = "right",
  dark,
  buttonDataTestId = "",
}: Props) => {
  return (
    <Tooltip
      sdsStyle={dark ? "dark" : "light"}
      placement={placement}
      width="wide"
      arrow
      title={<StyledTooltip>{text}</StyledTooltip>}
    >
      <TooltipButton
        data-testid={buttonDataTestId}
        sdsStyle="minimal"
        sdsType="secondary"
        isAllCaps={false}
      >
        <StyledIconImage src={questionMarkIcon} />
      </TooltipButton>
    </Tooltip>
  );
};

export default HelpTooltip;
