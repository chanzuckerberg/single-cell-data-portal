import styled from "@emotion/styled";
import { css } from "@emotion/react";
import { HEADER_HEIGHT_PX } from "src/components/Header/style";

export const TOP_PADDING_PX = 32;
export const CELL_GUIDE_SIDE_BAR_WIDTH_PX = 240;

export const CellGuideSidebarWrapper = styled.div`
  width: ${CELL_GUIDE_SIDE_BAR_WIDTH_PX}px;
`;

export const StickyWrapper = styled.div`
  display: block;
  position: sticky;
  top: ${HEADER_HEIGHT_PX + TOP_PADDING_PX}px;
`;

export const SearchBarWrapper = styled.div`
  margin-bottom: 64px;
  width: 240px;
`;

export const StickySidebarStyle = css`
  /* Hack because main has a global overflow CSS prop which interferes with sticky sidebar */
  main {
    overflow: unset !important;
  }
`;
