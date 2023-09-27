import { useCallback, useContext, useState } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { StateContext } from "src/views/WheresMyGene/common/store";
import { CITATION_NOTIFICATION_TEXT } from "src/views/WheresMyGeneV2/common/constants";

export const useConnect = () => {
  const state = useContext(StateContext);
  const { selectedGenes } = state;
  const [showCitationNotification, setShowCitationNotification] = useState(0);
  const copyCitation = useCallback(() => {
    track(EVENTS.WMG_CITATION_CLICKED);
    setShowCitationNotification((prev) => prev + 1);
    navigator.clipboard.writeText(CITATION_NOTIFICATION_TEXT);
  }, []);

  return {
    selectedGenes,
    showCitationNotification,
    copyCitation,
  };
};
