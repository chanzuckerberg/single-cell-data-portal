import { useCallback, useContext, useState } from "react";
import { StateContext } from "src/views/WheresMyGeneV2/common/store";
import { useFilterDimensions } from "src/common/queries/wheresMyGene";
import { EVENTS } from "src/common/analytics/events";
import { track } from "src/common/analytics";
import { COLOR_LEGEND } from "../RelativeGeneExpression/constants";

export const useConnect = () => {
  const [hoverStartTime, setHoverStartTime] = useState(0);
  const useHandleHoverEnd = (event: EVENTS, payload = {}) => {
    return useCallback(() => {
      if (Date.now() - hoverStartTime > 2 * 1000) {
        track(event, payload);
      }
    }, [event, payload]);
  };

  const handleLegendHoverEnd = useHandleHoverEnd(
    EVENTS.WMG_LEGEND_QUESTION_BUTTON_HOVER,
    { label: COLOR_LEGEND }
  );

  const { selectedFilters } = useContext(StateContext);
  const { data: filterDimensions } = useFilterDimensions();

  let { datasets = [] } = filterDimensions;

  if (selectedFilters.datasets.length > 0) {
    datasets = datasets.filter((dataset) =>
      selectedFilters.datasets.includes(dataset.id)
    );
  }

  const referenceCount = datasets.length;

  return { referenceCount, setHoverStartTime, handleLegendHoverEnd };
};
