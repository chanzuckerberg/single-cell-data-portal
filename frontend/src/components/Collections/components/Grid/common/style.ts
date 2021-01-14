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

export const CollectionHeaderCell = styled.th`
  ${titleColWidthCSS}
  text-align: left !important;
  padding: 0 !important;
`;

export const DatasetHeaderCell = styled.th`
  ${datasetTitleColWidthCSS}
  text-align: left !important;
  padding: 0 !important;
  padding-left: ${PT_GRID_SIZE_PX * 2}px !important;
`;

export const LeftAlignedHeaderCell = styled.th`
  ${detailsColWidthCSS}
  text-align: left !important;
  padding: 0 !important;
  padding-left: ${PT_GRID_SIZE_PX * 2}px !important;
`;

export const RightAlignedHeaderCell = styled.th`
  ${datasetTitleColWidthCSS}
  text-align: right !important;
  padding: 0 !important;
  padding-left: ${PT_GRID_SIZE_PX * 2}px !important;
`;
