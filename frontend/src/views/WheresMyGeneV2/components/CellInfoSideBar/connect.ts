import { useCallback, useContext, useState } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { useMarkerGenes } from "src/common/queries/wheresMyGene";
import { DispatchContext, State } from "../../../WheresMyGene/common/store";
import { addSelectedGenes } from "../../../WheresMyGene/common/store/actions";
import { MARKER_GENE_LABEL, MARKER_SCORE_LABEL } from "./constants";

export const useConnect = ({
  cellInfoCellType,
}: {
  cellInfoCellType: Exclude<State["cellInfoCellType"], null>;
}) => {
  const dispatch = useContext(DispatchContext);
  const [hoverStartTime, setHoverStartTime] = useState(0);

  const urlParams = new URLSearchParams(window.location.search);
  let testType: "binomtest" | undefined = undefined;
  if (urlParams.get("test") === "binomtest") {
    testType = "binomtest";
  }
  const { isLoading, data } = useMarkerGenes({
    cellTypeID: cellInfoCellType.cellType.id,
    organismID: cellInfoCellType.organismID,
    test: testType,
    tissueID: cellInfoCellType.tissueID,
  });

  const handleCopyGenes = useCallback(() => {
    if (!data) return;
    const genes = Object.keys(data.marker_genes);
    navigator.clipboard.writeText(genes.join(", "));
    track(EVENTS.WMG_FMG_COPY_GENES_CLICKED);
  }, [data]);

  const handleDisplayGenes = useCallback(() => {
    if (!data || !dispatch) return;
    const genes = Object.keys(data.marker_genes);
    dispatch(addSelectedGenes(genes));
    track(EVENTS.WMG_FMG_ADD_GENES_CLICKED);
  }, [data, dispatch]);

  const useHandleHoverEnd = (event: EVENTS, payload = {}) => {
    return useCallback(() => {
      if (Date.now() - hoverStartTime > 2 * 1000) {
        track(event, payload);
      }
    }, [event, payload]);
  };

  const handleFmgHoverEnd = useHandleHoverEnd(
    EVENTS.WMG_FMG_QUESTION_BUTTON_HOVER,
    { label: MARKER_GENE_LABEL }
  );
  const handleMarkerScoreHoverEnd = useHandleHoverEnd(
    EVENTS.WMG_FMG_QUESTION_BUTTON_HOVER,
    { label: MARKER_SCORE_LABEL }
  );

  return {
    handleCopyGenes,
    isLoading,
    handleDisplayGenes,
    data,
    handleFmgHoverEnd,
    handleMarkerScoreHoverEnd,
    setHoverStartTime,
  };
};
