import styled from "@emotion/styled";
import { Notification } from "@czi-sds/components";
import { NotificationWrapper } from "src/components/common/Filter/common/style";
import { beta100, beta400, beta600 } from "src/common/theme";

export const BETA_NOTIFICATION_ID = "beta-notification";

export const StyledNotification = styled(Notification)`
  /*
   * beta intent does not exist for SDS Notification, but the colors do targeting
   * specific id because using the generic .elevated class selector was messing with the copy link notification
   **/
  #beta-notification {
    border-color: ${beta400} !important;
    background-color: ${beta100};

    path {
      fill: ${beta600};
    }
  }
`;

export const StyledNotificationWrapper = styled(NotificationWrapper)`
  bottom: 10px;
  right: 30px;
  width: 360px;
`;

export const SubmitIssue = styled.a`
  color: ${beta600};
`;
