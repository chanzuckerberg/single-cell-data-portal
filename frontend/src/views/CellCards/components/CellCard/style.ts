import styled from "@emotion/styled";
import { fontHeaderXxl, Tag } from "czifui";

export const TOP_PADDING_PX = 32;

export const CellCardsView = styled.div`
  display: flex;
  flex-direction: row;
  column-gap: 120px;
  margin: auto;
  max-width: 1440px;
  padding: ${TOP_PADDING_PX}px 40px 40px 40px;
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
