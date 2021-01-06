import { HTMLTable } from "@blueprintjs/core";
import { GRAY } from "src/components/common/theme";
import styled, { css } from "styled-components";

export const textClippingCSS = css`
  white-space: nowrap;
`;

export const detailsColWidthCSS = css`
  width: calc(1 / 8 * 100%);
`;

export const StyledCell = styled.td`
  vertical-align: middle !important;
`;

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

// ROW

export const DetailsCell = styled(StyledCell)`
  color: ${GRAY.A};
  ${textClippingCSS}
  ${detailsColWidthCSS}
`;

export const LeftAlignedDetailsCell = styled(DetailsCell)`
  text-align: left !important;
`;

export const StyledRow = styled.tr`
  align-content: center;
  vertical-align: middle;
`;

export const RightAlignedDetailsCell = styled(DetailsCell)`
  text-align: right !important;
`;
