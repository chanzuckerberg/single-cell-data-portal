import { Notification } from "@czi-sds/components";
import {
  StyledNotificationWrapper,
  StyledNotificationLabel,
  StyledIcon,
  StyledNotificationDetails,
  StyledButtonDiv,
  StyledLabel,
  StyledButtonIcon,
} from "./style";
import { noop } from "src/common/constants/utils";
import { EXCLUDE_IN_SCREENSHOT_CLASS_NAME } from "../SaveExport";
import { useConnect } from "./connect";
import {
  CITATION_BUTTON_LABEL,
  CITATION_NOTIFICATION_LABEL,
  CITATION_NOTIFICATION_TEXT,
} from "src/views/WheresMyGeneV2/common/constants";

export default function CitationButton(): JSX.Element {
  const { copyCitation, showCitationNotification, selectedGenes } =
    useConnect();

  return (
    <>
      {showCitationNotification > 0 && (
        <StyledNotificationWrapper>
          <Notification
            key={showCitationNotification}
            autoDismiss={5000}
            onClose={noop}
            slideDirection="left"
            intent="info"
            icon={
              <StyledIcon sdsIcon={"copy"} sdsSize={"s"} sdsType={"static"} />
            }
          >
            <StyledNotificationLabel data-testid="citation-link-notification">
              {CITATION_NOTIFICATION_LABEL}
            </StyledNotificationLabel>
            <StyledNotificationDetails>
              {CITATION_NOTIFICATION_TEXT}
            </StyledNotificationDetails>
          </Notification>
        </StyledNotificationWrapper>
      )}
      <StyledButtonDiv className={EXCLUDE_IN_SCREENSHOT_CLASS_NAME}>
        <StyledLabel>{CITATION_BUTTON_LABEL}</StyledLabel>
        <StyledButtonIcon
          data-testid={"citation-button"}
          onClick={copyCitation}
          sdsSize="medium"
          sdsType="primary"
          sdsIcon="quote"
          disabled={selectedGenes.length === 0}
        />
      </StyledButtonDiv>
    </>
  );
}
