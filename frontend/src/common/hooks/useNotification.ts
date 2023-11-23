import { useContext } from "react";
import uuid from "react-uuid";
import { DispatchContext } from "src/views/WheresMyGeneV2/common/store";
import {
  addNotification,
  //   clearNotification,
} from "src/views/WheresMyGeneV2/common/store/actions";
import { ExposedNotificationProps } from "src/views/WheresMyGeneV2/common/store/reducer";

export function useNotification() {
  const dispatch = useContext(DispatchContext);

  function createNotification({
    message,
    intent,
    sdsIcon,
    label,
  }: ExposedNotificationProps) {
    const notificationId = uuid();

    if (!dispatch) return;
    dispatch(
      addNotification({ message, notificationId, intent, sdsIcon, label })
    );
    // setTimeout(() => {
    //   dispatch(clearNotification({ notificationId }));
    // }, 10000);
  }

  return {
    createNotification,
  };
}
