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
    <StyledButtonDiv >
      <StyledLabel>Source Data</StyledLabel>

      <StyledIconButton
        data-test-id={"source-data-button"}
        onClick={handleRightSidebarButtonClick}
        {...{
          // (thuang): Move this back to explicit prop={value} after
          // upgrading SDS to enable type checking again
          disabled: false,
          sdsSize: "medium",
          sdsType: "primary",
        }}
      >
        <Icon sdsIcon="infoCircle" sdsSize="l" sdsType="iconButton" />
      </StyledIconButton>
    </StyledButtonDiv>
  );
}
