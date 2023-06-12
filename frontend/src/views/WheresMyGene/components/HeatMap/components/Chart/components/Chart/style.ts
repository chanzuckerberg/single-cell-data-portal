import styled from "@emotion/styled";
import { MARGIN_BETWEEN_HEATMAPS } from "src/views/WheresMyGene/common/constants";

export const ChartContainer = styled.div`
  ${getWidthAndHeight}
  margin-bottom: ${MARGIN_BETWEEN_HEATMAPS}px;
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
