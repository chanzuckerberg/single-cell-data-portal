import styled from "@emotion/styled";
import { X_AXIS_CHART_HEIGHT_PX } from "src/views/WheresMyGene/common/constants";

export const ChartContainer = styled.div`
  ${getWidthAndHeight}
  margin-bottom: ${X_AXIS_CHART_HEIGHT_PX}px;
`;

function getWidthAndHeight({
  width,
  height,
}: {
  width: number;
  height: number;
}) {
  return `
    width: ${width}px;
    height: ${height}px;
  `;
}
