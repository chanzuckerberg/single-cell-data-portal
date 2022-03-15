import { css } from "@emotion/css";
/* eslint-disable @typescript-eslint/ban-ts-comment --- TODO revisit lint error below */
// @ts-ignore - TODO revisit lint error
import { TooltipTable } from "czifui";
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

export const StyledTooltipTable = styled(TooltipTable)``;

export const tooltipCss = css`
  margin: 0;
  margin-top: 20px;
  max-width: 400px !important;
`;
