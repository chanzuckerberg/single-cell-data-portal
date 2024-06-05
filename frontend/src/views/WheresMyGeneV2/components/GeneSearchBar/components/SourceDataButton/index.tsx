import { MouseEventHandler } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { EXCLUDE_IN_SCREENSHOT_CLASS_NAME } from "../SaveExport";
import {
  BadgeCounter,
  StyledButtonDiv,
  StyledButtonIcon,
  StyledLabel,
} from "./style";

export default function SourceDataButton({
  handleRightSidebarButtonClick,
  datasetsCount,
}: {
  handleRightSidebarButtonClick: MouseEventHandler<HTMLButtonElement>;
  datasetsCount: number;
}): JSX.Element {
  return (
    <StyledButtonDiv className={EXCLUDE_IN_SCREENSHOT_CLASS_NAME}>
      <StyledLabel>Source Data</StyledLabel>
      {datasetsCount > 0 && (
        <BadgeCounter
          badgeContent={datasetsCount}
          width={datasetsCount > 99 ? "28px" : "20px"}
          data-testid="source-data-badge"
        />
      )}
      <StyledButtonIcon
        data-testid={"source-data-button"}
        onClick={(event) => {
          track(EVENTS.WMG_SOURCE_DATA_CLICKED);
          handleRightSidebarButtonClick(event);
        }}
        sdsStyle="icon"
        sdsSize="medium"
        sdsType="primary"
        icon="InfoCircle"
      />
    </StyledButtonDiv>
  );
}
