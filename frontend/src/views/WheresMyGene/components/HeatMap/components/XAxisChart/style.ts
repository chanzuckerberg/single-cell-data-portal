import styled from "styled-components";
import { X_AXIS_CHART_HEIGHT_PX, Y_AXIS_CHART_WIDTH_PX } from "../../utils";

export const XAxisWrapper = styled.div`
  ${xAxisWidth}

  background-color: white;
  height: ${X_AXIS_CHART_HEIGHT_PX}px;
  position: sticky;
  top: 0;
  z-index: 2;
`;

export const XAxisMask = styled.div`
  background-color: white;
  height: ${X_AXIS_CHART_HEIGHT_PX}px;
  width: ${Y_AXIS_CHART_WIDTH_PX - 5}px;
  position: sticky;
  left: 0;
`;

export const XAxisContainer = styled.div`
  ${xAxisWidth}
  z-index: -1;
  background-color: white;
  height: ${X_AXIS_CHART_HEIGHT_PX}px;
  position: absolute;
`;
export const CellCountLabel = styled.div`
  font: 12px sans-serif;
  color: #6e7079;
  width: 12px;
  text-orientation: sideways;
  writing-mode: vertical-rl;
  position: sticky;
  left: 280px;
  padding-top: 16px;
  z-index: 2;
`;
function xAxisWidth({ width }: { width: number }) {
  return `
    width: ${width}px;
  `;
}
