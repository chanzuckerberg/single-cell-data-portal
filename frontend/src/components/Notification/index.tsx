import { Notification as SDSNotification } from "@czi-sds/components";
import { noop } from "src/common/constants/utils";
import {
  StyledIcon,
  StyledNotificationDetails,
  StyledNotificationLabel,
} from "./style";
import { NotificationWrapper } from "../common/Filter/common/style";
import { useContext } from "react";
import { StateContext } from "src/views/WheresMyGeneV2/common/store";
import { NOTIFICATION_AUTO_DISMISS_TIMEOUT } from "./constants";

export default function Notification(): JSX.Element {
  const state = useContext(StateContext);
  const { notifications } = state;
  return (
    <NotificationWrapper>
      {notifications.length > 0 &&
        notifications.map((notification) => (
          <SDSNotification
            key={notification.notificationId}
            autoDismiss={NOTIFICATION_AUTO_DISMISS_TIMEOUT}
            onClose={noop}
            slideDirection="left"
            intent={notification.intent}
            icon={
              <StyledIcon
                sdsIcon={notification.sdsIcon}
                sdsSize={notification.sdsSize}
                sdsType={"static"}
              />
            }
          >
            <StyledNotificationLabel data-testid="notification">
              {notification.label}
            </StyledNotificationLabel>
            <StyledNotificationDetails
              isCitation={notification?.isCitation ?? false}
            >
              {notification.message}
            </StyledNotificationDetails>
          </SDSNotification>
        ))}
    </NotificationWrapper>
  );
}
