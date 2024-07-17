import { useCallback, useState } from "react";
import {
  Notification,
  UseCollectionNotification,
} from "src/views/Collection/hooks/useNotification/types";
import uuid from "react-uuid";
import { NOTIFICATION_PROPS } from "src/views/Collection/hooks/useNotification/constants";

/**
 * Notification functionality for collection view notification component.
 * @returns notification functionality.
 */
export function useCollectionNotification(): UseCollectionNotification {
  const [notification, setNotification] = useState<Notification>();

  const onNotify = useCallback(
    (nextNotification: Omit<Notification, "id">): void => {
      setNotification({
        ...NOTIFICATION_PROPS,
        ...nextNotification,
        id: uuid(), // Use a unique id for each notification to ensure any new notifications created are displayed.
      });
    },
    []
  );

  return {
    notification,
    onNotify,
  };
}
