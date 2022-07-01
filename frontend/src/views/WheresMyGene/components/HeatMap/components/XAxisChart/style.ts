import styled from "styled-components";
import { X_AXIS_CHART_HEIGHT_PX, Y_AXIS_CHART_WIDTH_PX } from "../../utils";

export const XAxisWrapper = styled.div`
  ${xAxisWidth}

  background-color: white;
  height: ${X_AXIS_CHART_HEIGHT_PX}px;
  position: relative;
  top: 0;
  z-index: 2;
`;

export const XAxisMask = styled.div`
  background-color: white;
  height: ${X_AXIS_CHART_HEIGHT_PX}px;
  width: ${Y_AXIS_CHART_WIDTH_PX}px;
  position: relative;
`;

export const XAxisContainer = styled.div`
  ${xAxisWidth}

  background-color: white;
  height: ${X_AXIS_CHART_HEIGHT_PX}px;
  position: relative;
`;

function xAxisWidth({ width }: { width: number }) {
  return `
    width: ${width}px;
  `;
}
