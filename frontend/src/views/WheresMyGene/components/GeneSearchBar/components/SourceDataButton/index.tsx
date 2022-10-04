import { Icon } from "czifui";
import { MouseEventHandler } from "react";
import { StyledIconButton } from "../QuickSelect/style";
import { StyledButtonDiv, StyledLabel } from "./style";

export default function SourceDataButton({
  handleRightSidebarButtonClick,
}: {
  handleRightSidebarButtonClick: MouseEventHandler<HTMLButtonElement>;
}): JSX.Element {
  return (
    <StyledButtonDiv className="screenshot-exclude">
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
