import { useCallback, useState } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import {
  CITATION_NOTIFICATION_LABEL,
  CITATION_NOTIFICATION_TEXT,
} from "./constants";
import { useNotification } from "src/common/hooks/useNotification";
import { NOTIFICATION_AUTO_DISMISS_TIMEOUT } from "src/components/Notification/constants";

export const useConnect = () => {
  const { createNotification } = useNotification();
  const [isDisabled, setIsDisabled] = useState(false);
  const copyCitation = useCallback(() => {
    setIsDisabled(true);
    track(EVENTS.WMG_CITATION_CLICKED);
    createNotification({
      message: CITATION_NOTIFICATION_TEXT,
      intent: "info",
      sdsIcon: "quote",
      sdsSize: "l",
      label: CITATION_NOTIFICATION_LABEL,
      isCitation: true,
    });
    navigator.clipboard.writeText(CITATION_NOTIFICATION_TEXT);
    setTimeout(() => {
      setIsDisabled(false);
    }, NOTIFICATION_AUTO_DISMISS_TIMEOUT);
  }, [createNotification]);

  return {
    copyCitation,
    isDisabled,
  };
};
