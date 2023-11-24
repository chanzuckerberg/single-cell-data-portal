import { useContext } from "react";
import uuid from "react-uuid";
import { DispatchContext } from "src/views/WheresMyGeneV2/common/store";
import {
  addNotification,
  //   clearNotification,
} from "src/views/WheresMyGeneV2/common/store/actions";
import { NotificationProps } from "src/views/WheresMyGeneV2/common/store/reducer";

export function useNotification() {
  const dispatch = useContext(DispatchContext);

  function createNotification({
    message,
    intent,
    sdsIcon,
    label,
    sdsSize,
    isCitation,
  }: NotificationProps) {
    const notificationId = uuid();

    if (!dispatch) return;
    dispatch(
      addNotification({
        message,
        notificationId,
        intent,
        sdsIcon,
        sdsSize,
        label,
        isCitation,
      })
    );
  }

  return {
    createNotification,
  };
}
