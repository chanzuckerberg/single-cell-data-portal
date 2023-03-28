import { noop } from "src/common/constants/utils";
import { EXCLUDE_IN_SCREENSHOT_CLASS_NAME } from "../GeneSearchBar/components/SaveExport";
import {
  StyledNotification,
  StyledNotificationWrapper,
  SubmitIssue,
} from "./style";

export default function Beta(): JSX.Element {
  return (
    <StyledNotificationWrapper>
      <StyledNotification
        id="beta-notification"
        intent="info"
        autoDismiss={false}
        dismissDirection="left"
        onClose={noop}
        className={EXCLUDE_IN_SCREENSHOT_CLASS_NAME}
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
