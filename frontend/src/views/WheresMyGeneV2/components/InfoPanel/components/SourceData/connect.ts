import { useCallback, useContext, useMemo, useState } from "react";
import { StateContext } from "src/views/WheresMyGeneV2/common/store";
import {
  aggregateCollectionsFromDatasets,
  useFilterDimensions,
} from "src/common/queries/wheresMyGene";
import { Collections } from "./types";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { SOURCE_DATA } from "./constants";

export const useConnect = () => {
  const [hoverStartTime, setHoverStartTime] = useState(0);
  const useHandleHoverEnd = (event: EVENTS, payload = {}) => {
    return useCallback(() => {
      if (Date.now() - hoverStartTime > 2 * 1000) {
        track(event, payload);
      }
    }, [event, payload]);
  };

  const handleBadgeHoverEnd = useHandleHoverEnd(
    EVENTS.WMG_LEGEND_QUESTION_BUTTON_HOVER,
    { label: SOURCE_DATA }
  );

  const { selectedFilters } = useContext(StateContext);
  const { data: filterDimensions } = useFilterDimensions();

  let { datasets = [] } = filterDimensions;
  if (selectedFilters.datasets.length > 0) {
    datasets = datasets.filter((dataset) =>
      selectedFilters.datasets.includes(dataset.id)
    );
  }
  const collections: Collections = useMemo(() => {
    return aggregateCollectionsFromDatasets(datasets);
  }, [datasets]);

  return { collections, setHoverStartTime, handleBadgeHoverEnd };
};
