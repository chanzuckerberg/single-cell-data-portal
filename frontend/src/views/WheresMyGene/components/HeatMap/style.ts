import styled from "@emotion/styled";
import { X_AXIS_CHART_TOTAL_HEIGHT_PX } from "./components/XAxisChart/style";
import { Y_AXIS_CHART_WIDTH_PX } from "./utils";

export const CHART_LEFT_PADDING = 10;

export const Container = styled.div`
  height: 75vh;
  width: 100%;
  overflow: scroll;
  position: relative;
`;

export const YAxisWrapper = styled.div`
  width: ${Y_AXIS_CHART_WIDTH_PX}px;
  position: absolute;
  top: ${X_AXIS_CHART_TOTAL_HEIGHT_PX}px;
  left: 0;
  z-index: 1;
  padding-top: 5px;
  /* Somehow Firefox requires this to scroll */
  overflow: hidden;

  ${({ height }: { height: number }) => {
    return `
      height: ${height}px
    `;
  }}
`;

export const XAxisMask = styled.div`
  width: ${Y_AXIS_CHART_WIDTH_PX}px;
  height: ${X_AXIS_CHART_TOTAL_HEIGHT_PX}px;
`;
export const XAxisWrapper = styled.div`
  display: flex;
  background-color: white;
  flex-direction: row;
  position: sticky;
  top: 0;
  width: 100%;
  z-index: 2;
`;

export const ChartRowWrapper = styled.div`
  display: flex;
  flex-direction: row;
  position: relative;
  left: ${Y_AXIS_CHART_WIDTH_PX}px;
  column-gap: 40px;
`;

export const ChartWrapper = styled.div`
  position: absolute;
  padding-left: ${CHART_LEFT_PADDING}px;
  padding-top: 5px;
  top: ${X_AXIS_CHART_TOTAL_HEIGHT_PX}px;
`;
