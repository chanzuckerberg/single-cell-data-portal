import {
  ExposedNotificationProps,
  NotificationProps,
} from "@czi-sds/components";

export type CollectionNotification = UseCollectionNotification;

export interface Notification extends NotificationProps {
  id: string;
  onClose?: ExposedNotificationProps["onClose"];
  title: string;
}

export type OnNotifyFn = (notification: Omit<Notification, "id">) => void;

export interface UseCollectionNotification {
  notification: Notification | undefined;
  onNotify: OnNotifyFn;
}
