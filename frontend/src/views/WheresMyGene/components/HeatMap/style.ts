import styled from "styled-components";
import { X_AXIS_CHART_HEIGHT, Y_AXIS_CHART_WIDTH } from "./utils";
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

export const ChartContainer = styled.div`
  left: 0px;
  top: 0px;
  position: absolute;

  ${({ width, height }: { width: number; height: number }) => {
    return `
      width: ${width}px;
      height: ${height}px;
    `;
  }})}
`;

export const XAxisWrapper = styled.div`
  background-color: white;
  height: ${X_AXIS_CHART_HEIGHT};
  position: sticky;
  top: 0;
  z-index: 2;

  ${xAxisWidth}
`;

export const XAxisMask = styled.div`
  background-color: white;
  height: ${X_AXIS_CHART_HEIGHT};
  width: ${Y_AXIS_CHART_WIDTH};
  position: sticky;
  left: 0;
`;

export const XAxisContainer = styled.div`
  background-color: white;
  height: ${X_AXIS_CHART_HEIGHT};
  position: absolute;

  ${xAxisWidth}
`;

export const YAxisContainer = styled.div`
  background-color: white;
  width: ${Y_AXIS_CHART_WIDTH};
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

function xAxisWidth({ width }: { width: number }) {
  return `
    width: ${width}px;
  `;
}
