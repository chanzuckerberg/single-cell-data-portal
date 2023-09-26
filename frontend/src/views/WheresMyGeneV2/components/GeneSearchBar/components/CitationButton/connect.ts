import { useCallback, useContext, useState } from "react";
import { StateContext } from "src/views/WheresMyGene/common/store";
import { CITATION_TEXT } from "src/views/WheresMyGeneV2/common/constants";

export const useConnect = () => {
  const state = useContext(StateContext);
  const { selectedGenes } = state;
  const [showCitationNotification, setShowCitationNotification] = useState(0);
  const copyCitation = useCallback(() => {
    setShowCitationNotification((prev) => prev + 1);
    navigator.clipboard.writeText(CITATION_TEXT);
  }, []);

  return {
    selectedGenes,
    showCitationNotification,
    copyCitation,
  };
};
