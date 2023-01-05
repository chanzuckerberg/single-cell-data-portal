import { MouseEventHandler } from "react";
import { StyledButtonIcon } from "../QuickSelect/style";
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

      <StyledButtonIcon
        data-test-id={"source-data-button"}
        onClick={handleRightSidebarButtonClick}
        sdsSize="medium"
        sdsType="primary"
        sdsIcon="infoCircle"
      />
    </StyledButtonDiv>
  );
}
