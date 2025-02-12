import { useCallback, useContext, useMemo, useState } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { useMarkerGenes } from "src/common/queries/wheresMyGene";
import { DispatchContext, State, StateContext } from "../../common/store";
import { addSelectedGenes } from "../../common/store/actions";
import { HOVER_START_TIME_MS } from "../../common/constants";
import {
  MARKER_GENE_LABEL,
  MARKER_SCORE_LABEL,
  SPECIFICITY_LABEL,
} from "src/common/constants/markerGenes";
import { generateDifferentialExpressionUrl } from "../GeneSearchBar/components/ShareButton/utils";

export const useConnect = ({
  cellInfoCellType,
}: {
  cellInfoCellType: Exclude<State["cellInfoCellType"], null>;
}) => {
  const { selectedFilters, selectedOrganismId } = useContext(StateContext);
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
      if (Date.now() - hoverStartTime > HOVER_START_TIME_MS) {
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
  const handleSpecificityHoverEnd = useHandleHoverEnd(
    EVENTS.WMG_FMG_QUESTION_BUTTON_HOVER,
    { label: SPECIFICITY_LABEL }
  );

  const differentialExpressionUrl = useMemo(
    () =>
      generateDifferentialExpressionUrl({
        filters: selectedFilters,
        organism: selectedOrganismId,
        tissue: cellInfoCellType.tissueID,
        cellType: cellInfoCellType.cellType.id,
      }),
    [
      selectedFilters,
      selectedOrganismId,
      cellInfoCellType.tissueID,
      cellInfoCellType.cellType.id,
    ]
  );

  return {
    handleCopyGenes,
    isLoading,
    handleDisplayGenes,
    data,
    handleFmgHoverEnd,
    handleMarkerScoreHoverEnd,
    handleSpecificityHoverEnd,
    setHoverStartTime,
    differentialExpressionUrl,
  };
};
