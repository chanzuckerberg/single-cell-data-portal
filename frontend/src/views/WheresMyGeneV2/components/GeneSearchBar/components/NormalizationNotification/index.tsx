import { Notification } from "@czi-sds/components";
import { noop } from "src/common/constants/utils";
import {
  StyledIcon,
  StyledLink,
  StyledLinkBody,
  StyledLinkWrapper,
  StyledNotificationDetails,
  StyledNotificationLabel,
  StyledNotificationWrapper,
} from "../ShareButton/style";
import { useConnect } from "./connect";
import {
  NOMALIZATION_NOTIFICATION_BODY_LINK,
  NORMALIZATION_NOTIFICATION_BODY,
  NORMALIZATION_NOTIFICATION_BODY_2,
  NORMALIZATION_NOTIFICATION_BODY_LINK_LABEL,
  NORMALIZATION_NOTIFICATION_DOCS_URL,
  NORMALIZATION_NOTIFICATION_LINK_TEXT,
  NORMALIZATION_NOTIFICATION_TITLE,
} from "./constants";

export default function NormalizationNotification(): JSX.Element {
  const { showURLCopyNotification } = useConnect();

  return (
    <StyledNotificationWrapper>
      <Notification
        key={showURLCopyNotification}
        autoDismiss={10000}
        onClose={noop}
        slideDirection={"left"}
        intent={"info"}
        icon={<StyledIcon sdsIcon={"link"} sdsSize={"s"} sdsType={"static"} />}
      >
        <StyledNotificationLabel data-testid="share-link-notification">
          {NORMALIZATION_NOTIFICATION_TITLE}
        </StyledNotificationLabel>
        <StyledNotificationDetails>
          {NORMALIZATION_NOTIFICATION_BODY}
          <StyledLinkBody href={NOMALIZATION_NOTIFICATION_BODY_LINK}>
            {NORMALIZATION_NOTIFICATION_BODY_LINK_LABEL}
          </StyledLinkBody>
          {NORMALIZATION_NOTIFICATION_BODY_2}
          <StyledLinkWrapper>
            <StyledLink href={NORMALIZATION_NOTIFICATION_DOCS_URL}>
              {NORMALIZATION_NOTIFICATION_LINK_TEXT}
            </StyledLink>
          </StyledLinkWrapper>
        </StyledNotificationDetails>
      </Notification>
    </StyledNotificationWrapper>
  );
}
