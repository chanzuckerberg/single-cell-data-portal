import styled from "styled-components";

export const ChartContainer = styled.div`
  ${getWidthAndHeight}
`;

export const Wrapper = styled.div`
  ${getWidthAndHeight}
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
