import styled, { css } from "styled-components";

export const textClippingCSS = css`
  white-space: nowrap;
  text-overflow: clip;
  overflow: hidden;
`;

export const detailsColWidthCSS = css`
  width: calc(1 / 8 * 100%);
`;

export const StyledCell = styled.td`
  vertical-align: middle !important;
`;
