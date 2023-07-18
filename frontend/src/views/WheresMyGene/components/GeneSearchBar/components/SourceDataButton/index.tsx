import { MouseEventHandler } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { StyledButtonIcon } from "../QuickSelect/style";
import { StyledButtonDiv, StyledLabel } from "./style";
import { EXCLUDE_IN_SCREENSHOT_CLASS_NAME } from "src/views/WheresMyGeneV2/components/GeneSearchBar/components/SaveExport";

export default function SourceDataButton({
  handleRightSidebarButtonClick,
}: {
  handleRightSidebarButtonClick: MouseEventHandler<HTMLButtonElement>;
}): JSX.Element {
  return (
    <StyledButtonDiv className={EXCLUDE_IN_SCREENSHOT_CLASS_NAME}>
      <StyledLabel>Source Data</StyledLabel>

      <StyledButtonIcon
        data-testid={"source-data-button"}
        onClick={(event) => {
          track(EVENTS.WMG_SOURCE_DATA_CLICKED);
          handleRightSidebarButtonClick(event);
        }}
        sdsSize="medium"
        sdsType="primary"
        sdsIcon="infoCircle"
      />
    </StyledButtonDiv>
  );
}
