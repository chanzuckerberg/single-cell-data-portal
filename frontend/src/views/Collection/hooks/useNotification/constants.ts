import { Notification } from "src/views/Collection/hooks/useNotification/types";
import { noop } from "src/common/constants/utils";

export const NOTIFICATION_PROPS: Partial<Notification> = {
  autoDismiss: true,
  onClose: noop,
};
