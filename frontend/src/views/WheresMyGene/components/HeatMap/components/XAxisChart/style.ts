import styled from "styled-components";
import { X_AXIS_CHART_HEIGHT_PX, Y_AXIS_CHART_WIDTH_PX } from "../../utils";

export const XAxisWrapper = styled.div`
  ${xAxisWidth}

  background-color: white;
  height: ${X_AXIS_CHART_HEIGHT_PX}px;
  display: inline-block;
  position: sticky;
  top: 0px;
  z-index: 2;
`;

export const XAxisMask = styled.div`
  background-color: white;
  height: ${X_AXIS_CHART_HEIGHT_PX}px;
  width: ${Y_AXIS_CHART_WIDTH_PX}px;
  position: sticky;
  display: inline-block;
  z-index: 2;
`;

export const XAxisContainer = styled.div`
  ${xAxisWidth}

  background-color: white;
  height: ${X_AXIS_CHART_HEIGHT_PX}px;
  position: sticky;
`;

function xAxisWidth({ width }: { width: number }) {
  return `
    width: ${width}px;
  `;
}
