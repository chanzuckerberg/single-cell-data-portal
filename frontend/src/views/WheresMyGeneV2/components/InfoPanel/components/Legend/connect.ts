import { useCallback, useState } from "react";
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

  const { data: filterDimensions } = useFilterDimensions();

  const { datasets = [] } = filterDimensions;

  const referenceCount = datasets.length;

  return { referenceCount, setHoverStartTime, handleLegendHoverEnd };
};
