import { HTMLTable } from "@blueprintjs/core";
import styled, { css } from "styled-components";
import { detailsColWidthCSS } from "./components/common/style";

const titleColWidthCSS = css`
  width: calc(3 / 8 * 100%);
`;

export const StyledCollectionsGrid = styled(HTMLTable)`
  grid-column: 1/9;
  margin-top: 16px;
`;

export const CollectionHeaderCell = styled.th`
  ${titleColWidthCSS}
  text-align: left !important;
`;

export const LeftAlignedHeaderCell = styled.th`
  ${detailsColWidthCSS}
  text-align: left !important;
  margin-left: 16px;
`;

export const RightAlignedHeaderCell = styled.th`
  width: calc(1 / 8 * 100%);
  text-align: right !important;
  margin-left: 16px;
`;
