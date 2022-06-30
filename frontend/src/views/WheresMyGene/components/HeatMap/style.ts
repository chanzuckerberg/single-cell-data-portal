import styled from "styled-components";
import { Y_AXIS_CHART_WIDTH_PX } from "./utils";

export const Container = styled.div`
  height: 75vh;
  width: 80vw;
  overflow: scroll;
  position: relative;
`;

export const YAxisWrapper = styled.div`
  width: ${Y_AXIS_CHART_WIDTH_PX}px;
  position: relative;
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
  position: relative;
  margin-top: -60px;
`;