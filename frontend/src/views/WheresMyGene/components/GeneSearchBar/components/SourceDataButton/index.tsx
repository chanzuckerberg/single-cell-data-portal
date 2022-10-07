import { Icon } from "czifui";
import { MouseEventHandler } from "react";
import { StyledIconButton } from "../QuickSelect/style";
import { EXCLUDE_IN_SCREENSHOT_CLASS_NAME } from "../SaveImage";
import { StyledButtonDiv, StyledLabel } from "./style";

export default function SourceDataButton({
  handleRightSidebarButtonClick,
}: {
  handleRightSidebarButtonClick: MouseEventHandler<HTMLButtonElement>;
}): JSX.Element {
  return (
    <StyledButtonDiv className={EXCLUDE_IN_SCREENSHOT_CLASS_NAME}>
      <StyledLabel>Source Data</StyledLabel>

      <StyledIconButton
        disabled={false}
        data-test-id={"download-button"}
        onClick={handleRightSidebarButtonClick}
        sdsType="primary"
        sdsSize="medium"
      >
        <Icon sdsIcon="infoCircle" sdsSize="l" sdsType="iconButton" />
      </StyledIconButton>
    </StyledButtonDiv>
  );
}
