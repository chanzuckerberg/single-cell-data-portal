import styled from "@emotion/styled";
import { css } from "@emotion/react";
import { HEADER_HEIGHT_PX } from "src/components/Header/style";
import { TOP_PADDING_PX } from "../../style";

export const CellGuideSidebarWrapper = styled.div`
  width: 240px;
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
