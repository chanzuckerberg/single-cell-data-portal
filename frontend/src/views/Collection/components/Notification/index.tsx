import { Props } from "src/views/Collection/components/Notification/types";
import { StyledNotificationWrapper } from "src/views/Collection/components/Notification/style";
import { Notification as SDSNotification } from "@czi-sds/components";

export default function Notification({
  notification,
}: Props): JSX.Element | null {
  if (!notification) return null;
  const { autoDismiss, id, intent, onClose, slideDirection, title } =
    notification || {};
  return (
    <StyledNotificationWrapper key={id}>
      <SDSNotification
        autoDismiss={autoDismiss}
        intent={intent}
        onClose={onClose}
        slideDirection={slideDirection}
      >
        {title}
      </SDSNotification>
    </StyledNotificationWrapper>
  );
}
