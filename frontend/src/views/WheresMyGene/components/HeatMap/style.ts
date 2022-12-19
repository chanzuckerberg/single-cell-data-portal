import styled from "@emotion/styled";
import { X_AXIS_CHART_HEIGHT_PX, Y_AXIS_CHART_WIDTH_PX } from "./utils";

export const Container = styled.div`
  height: 75vh;
  width: 100%;
  overflow: scroll;
  position: relative;
`;

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
  width: ${Y_AXIS_CHART_WIDTH_PX}px;
  position: sticky;
  left: 0;
`;

export const XAxisContainer = styled.div`
  ${xAxisWidth}

  background-color: white;
  height: ${X_AXIS_CHART_HEIGHT_PX}px;
  position: absolute;
`;

export const YAxisWrapper = styled.div`
  width: ${Y_AXIS_CHART_WIDTH_PX}px;
  position: sticky;
  top: 0;
  left: 0;
  z-index: 1;
  /* Somehow Firefox requires this to scroll */
  overflow: hidden;

  ${({ height }: { height: number }) => {
    return `
      height: ${height}px
    `;
  }}
`;

export const ChartWrapper = styled.div`
  display: flex;
  flex-direction: column;
  position: absolute;
  top: 0;
`;

function xAxisWidth({ width }: { width: number }) {
  return `
    width: ${width}px;
  `;
}
