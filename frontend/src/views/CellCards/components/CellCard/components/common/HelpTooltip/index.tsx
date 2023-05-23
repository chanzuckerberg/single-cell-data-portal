import { Tooltip } from "czifui";

import questionMarkIcon from "src/common/images/question-mark-icon.svg";
import { TooltipButton, StyledTooltip, StyledIconImage } from "./style";
import { ReactElement } from "react";

interface Props {
  text: string | ReactElement;
  placement?: "top" | "bottom" | "left" | "right";
}
const HelpTooltip = ({ text, placement = "right" }: Props) => {
  return (
    <Tooltip
      sdsStyle="light"
      placement={placement}
      width="wide"
      arrow
      title={<StyledTooltip>{text}</StyledTooltip>}
    >
      <TooltipButton sdsStyle="minimal" sdsType="secondary" isAllCaps={false}>
        <StyledIconImage src={questionMarkIcon} />
      </TooltipButton>
    </Tooltip>
  );
};

export default HelpTooltip;
