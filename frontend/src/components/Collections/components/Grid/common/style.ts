import { HTMLTable } from "@blueprintjs/core";
import { PT_GRID_SIZE_PX } from "src/components/common/theme";
import styled, { css } from "styled-components";

export const textClippingCSS = css`
  white-space: nowrap;
`;

export const detailsColWidthCSS = css`
  width: calc(1 / 8 * 100%);
`;

const titleColWidthCSS = css`
  width: calc(3 / 8 * 100%);
`;

const datasetTitleColWidthCSS = css`
  width: calc(2 / 8 * 100%);
`;

export const StyledCollectionsGrid = styled(HTMLTable)`
  grid-column: 1/9;
  margin-top: ${PT_GRID_SIZE_PX * 2}px;
`;

const HeaderCell = styled.th`
  padding: 12px 0 !important;
  font-weight: 500;
  font-size: 14px;
  line-height: 18px;
`;

export const CollectionHeaderCell = styled(HeaderCell)`
  ${titleColWidthCSS}
  text-align: left !important;
`;

export const DatasetHeaderCell = styled(HeaderCell)`
  ${datasetTitleColWidthCSS}
  text-align: left !important;
  padding-left: ${PT_GRID_SIZE_PX * 2}px !important;
`;

export const LeftAlignedHeaderCell = styled(HeaderCell)`
  ${detailsColWidthCSS}
  text-align: left !important;
  padding-left: ${PT_GRID_SIZE_PX * 2}px !important;
`;

export const RightAlignedHeaderCell = styled(HeaderCell)`
  ${detailsColWidthCSS}
  text-align: right !important;
  padding-left: ${PT_GRID_SIZE_PX * 2}px !important;
`;
