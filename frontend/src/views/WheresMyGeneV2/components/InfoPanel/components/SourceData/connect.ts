import { useCallback, useMemo, useState } from "react";
import {
  aggregateCollectionsFromDatasets,
  useFilterDimensions,
} from "src/common/queries/wheresMyGene";
import { Collections } from "./types";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { SOURCE_DATA } from "./constants";
import { EMPTY_ARRAY } from "src/common/constants/utils";
import { HOVER_START_TIME_MS } from "src/views/WheresMyGeneV2/common/constants";

export const useConnect = () => {
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

  const { data: filterDimensions } = useFilterDimensions();

  const { datasets = EMPTY_ARRAY } = filterDimensions;

  const collections: Collections = useMemo(() => {
    return aggregateCollectionsFromDatasets(datasets);
  }, [datasets]);

  return { collections, setHoverStartTime, handleBadgeHoverEnd };
};
