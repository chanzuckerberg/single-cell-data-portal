import { noop } from "src/common/constants/utils";

import {
  BETA_NOTIFICATION_ID,
  StyledNotification,
  StyledNotificationWrapper,
  SubmitIssue,
} from "./style";
import { EXCLUDE_IN_SCREENSHOT_CLASS_NAME } from "src/views/WheresMyGeneV2/components/GeneSearchBar/components/SaveExport";

export default function Beta(): JSX.Element {
  return (
    <StyledNotificationWrapper>
      <StyledNotification
        id={BETA_NOTIFICATION_ID}
        intent="info"
        autoDismiss={false}
        slideDirection="left"
        onClose={noop}
        className={EXCLUDE_IN_SCREENSHOT_CLASS_NAME}
        data-testid="survey-alert-id"
      >
        We would appreciate your feedback, please fill out a{" "}
        <SubmitIssue
          href="https://airtable.com/shrLwepDSEX1HI6bo"
          target="_blank"
          rel="noopener"
        >
          quick survey
        </SubmitIssue>
        .
      </StyledNotification>
    </StyledNotificationWrapper>
  );
}
