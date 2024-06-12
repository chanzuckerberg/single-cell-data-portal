import { useCallback, useContext, useMemo, useState } from "react";
import { aggregateCollectionsFromDatasets } from "src/common/queries/wheresMyGene";
import { Collections } from "./types";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { SOURCE_DATA } from "./constants";
import { EMPTY_ARRAY } from "src/common/constants/utils";
import { HOVER_START_TIME_MS } from "src/views/WheresMyGeneV2/common/constants";
import { StateContext } from "src/views/DifferentialExpression/common/store";
import { useQueryGroupFilterDimensions } from "src/common/queries/differentialExpression";
import { QueryGroups } from "src/views/DifferentialExpression/common/store/reducer";

export const useConnect = () => {
  const { submittedQueryGroups: queryGroupsRaw } = useContext(StateContext);

  // By construction, if SourceData is being rendered, then queryGroupsRaw is defined
  const queryGroups = queryGroupsRaw as QueryGroups;

  const [hoverStartTime, setHoverStartTime] = useState(0);
  const useHandleHoverEnd = (event: EVENTS, payload = {}) => {
    return useCallback(() => {
      if (Date.now() - hoverStartTime > HOVER_START_TIME_MS) {
        track(event, payload);
      }
    }, [event, payload]);
  };

  const handleBadgeHoverEnd = useHandleHoverEnd(
    EVENTS.WMG_LEGEND_QUESTION_BUTTON_HOVER,
    { label: SOURCE_DATA }
  );

  const { data: filterDimensions1 } = useQueryGroupFilterDimensions(
    queryGroups.queryGroup1
  );
  const { data: filterDimensions2 } = useQueryGroupFilterDimensions(
    queryGroups.queryGroup2
  );

  const { datasets: datasets1 = EMPTY_ARRAY } = filterDimensions1;
  const { datasets: datasets2 = EMPTY_ARRAY } = filterDimensions2;

  const collections1: Collections = useMemo(() => {
    return aggregateCollectionsFromDatasets(datasets1);
  }, [datasets1]);

  const collections2: Collections = useMemo(() => {
    return aggregateCollectionsFromDatasets(datasets2);
  }, [datasets2]);

  return { collections1, collections2, setHoverStartTime, handleBadgeHoverEnd };
};
