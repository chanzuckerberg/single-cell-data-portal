import { useCallback, useState } from "react";
import { useFilterDimensions } from "src/common/queries/wheresMyGene";
import { EVENTS } from "src/common/analytics/events";
import { track } from "src/common/analytics";
import { COLOR_LEGEND } from "../RelativeGeneExpression/constants";
import { HOVER_START_TIME } from "src/views/WheresMyGeneV2/common/constants";
import { EMPTY_ARRAY } from "src/common/constants/utils";

export const useConnect = () => {
  const [hoverStartTime, setHoverStartTime] = useState(0);
  const useHandleHoverEnd = (event: EVENTS, payload = {}) => {
    return useCallback(() => {
      if (Date.now() - hoverStartTime > HOVER_START_TIME) {
        track(event, payload);
      }
    }, [event, payload]);
  };

  const handleLegendHoverEnd = useHandleHoverEnd(
    EVENTS.WMG_LEGEND_QUESTION_BUTTON_HOVER,
    { label: COLOR_LEGEND }
  );

  const { data: filterDimensions } = useFilterDimensions();

  const { datasets = EMPTY_ARRAY } = filterDimensions;

  const referenceCount = datasets.length;

  return { referenceCount, setHoverStartTime, handleLegendHoverEnd };
};
