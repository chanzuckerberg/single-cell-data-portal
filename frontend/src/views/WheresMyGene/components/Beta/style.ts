import styled from "@emotion/styled";
import { CommonThemeProps, getColors, Notification } from "@czi-sds/components";
import { NotificationWrapper } from "src/components/common/Filter/common/style";

export const BETA_NOTIFICATION_ID = "beta-notification";

export const StyledNotification = styled(Notification)`
  ${(props: CommonThemeProps) => {
    const colors = getColors(props);

    // beta intent does not exist for SDS Notification, but the colors do
    // targeting specific id because using the generic .elevated class selector was messing with the copy link notification
    return `
      #beta-notification {
        border-color: ${colors?.beta[400]} !important;
        background-color: ${colors?.beta[100]};

        path {
          fill: ${colors?.beta[600]};
        }
      }
    `;
  }}
`;

export const StyledNotificationWrapper = styled(NotificationWrapper)`
  bottom: 10px;
  right: 30px;
  width: 360px;
`;

export const SubmitIssue = styled.a`
  ${(props: CommonThemeProps) => {
    const colors = getColors(props);

    return `
      color: ${colors?.beta[600]};
    `;
  }}
`;
