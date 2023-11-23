import { Notification as SDSNotification } from "@czi-sds/components";
import { noop } from "src/common/constants/utils";
import {
  StyledIcon,
  StyledNotificationDetails,
  StyledNotificationLabel,
} from "./style";
import { NotificationWrapper } from "../common/Filter/common/style";
import { Props } from "./types";
import { useContext } from "react";
import { StateContext } from "src/views/WheresMyGeneV2/common/store";

export default function Notification({}: Props): JSX.Element {
  const state = useContext(StateContext);
  const { notifications } = state;
  console.log(notifications);
  return (
    <NotificationWrapper>
      {notifications.length > 0 &&
        notifications.map((notification) => (
          <SDSNotification
            key={notification.notificationId}
            autoDismiss={5000}
            onClose={noop}
            slideDirection="left"
            intent={notification.intent}
            icon={
              <StyledIcon
                sdsIcon={notification.sdsIcon}
                sdsSize={"l"}
                sdsType={"static"}
              />
            }
          >
            <StyledNotificationLabel data-testid="notification">
              {notification.label}
            </StyledNotificationLabel>
            <StyledNotificationDetails>
              {notification.message}
            </StyledNotificationDetails>
          </SDSNotification>
        ))}
    </NotificationWrapper>
  );
}
