import { useCallback, useContext } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { StateContext } from "src/views/WheresMyGeneV2/common/store";
import {
  CITATION_NOTIFICATION_LABEL,
  CITATION_NOTIFICATION_TEXT,
} from "./constants";
import { useNotification } from "src/common/hooks/useNotification";

export const useConnect = () => {
  const { createNotification } = useNotification();
  const state = useContext(StateContext);
  const { selectedGenes } = state;
  const copyCitation = useCallback(() => {
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
  }, []);

  return {
    selectedGenes,
    copyCitation,
  };
};
