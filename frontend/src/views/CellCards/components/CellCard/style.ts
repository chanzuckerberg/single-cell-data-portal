import styled from "@emotion/styled";
import { fontHeaderXxl, Tag } from "@czi-sds/components";

export const TOP_PADDING_PX = 32;
export const SIDEBAR_COLUMN_GAP_PX = 120;
export const LEFT_RIGHT_PADDING_PX = 40;

export const CellCardsView = styled.div`
  display: flex;
  flex-direction: row;
  column-gap: ${SIDEBAR_COLUMN_GAP_PX}px;
  margin: auto;
  max-width: 1440px;
  padding: ${TOP_PADDING_PX}px ${LEFT_RIGHT_PADDING_PX}px 0px
    ${LEFT_RIGHT_PADDING_PX}px;
`;

export const Wrapper = styled.div`
  display: flex;
  flex-direction: column;
  align-self: stretch;
  min-width: 640px;
`;

export const CellCardHeaderInnerWrapper = styled.div`
  display: flex;
  column-gap: 8px;
  align-items: center;
`;

export const CellCardHeader = styled.div`
  display: flex;
  column-gap: 8px;
  flex-direction: row;
  align-items: center;
  justify-content: space-between;
`;

export const CellCardName = styled.div`
  ${fontHeaderXxl}
  font-weight: 700;
`;

export const StyledTag = styled(Tag)`
  height: 24px;
  margin: 0;
`;

export const SearchBarWrapper = styled.div`
  margin-bottom: 20px;
  width: 240px;
`;
