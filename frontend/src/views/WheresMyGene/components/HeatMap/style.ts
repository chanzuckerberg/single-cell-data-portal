import styled from "styled-components";
import { X_AXIS_CHART_HEIGHT_PX, Y_AXIS_CHART_WIDTH_PX } from "./utils";

export const Container = styled.div`
  height: 75vh;
  width: 80vw;
  overflow: scroll;
  position: relative;
`;

export const Loader = styled.div`
  position: fixed;
  top: 75px;
  left: 50vw;
  width: 200px;
  height: 50px;
  background-color: white;
  display: flex;
  align-items: center;
  gap: 10px;
  justify-content: center;
  box-shadow: 0px 0px 0px rgba(16, 22, 26, 0.1),
    0px 4px 8px rgba(16, 22, 26, 0.2);
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

  ${({ height }: { height: number }) => {
    return `
      height: ${height}px;
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
