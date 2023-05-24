import styled from "@emotion/styled";
import { HEADER_HEIGHT_PX } from "src/components/Header/style";
import { TOP_PADDING_PX } from "../../style";

export const CellCardsSidebarWrapper = styled.div`
  width: 240px;
`;

export const StickyWrapper = styled.div`
  display: block;
  position: sticky;
  top: ${HEADER_HEIGHT_PX + TOP_PADDING_PX}px;
`;

export const SearchBarWrapper = styled.div`
  margin-bottom: 20px;
  width: 240px;
`;
