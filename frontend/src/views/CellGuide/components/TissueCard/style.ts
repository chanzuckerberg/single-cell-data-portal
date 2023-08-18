import styled from "@emotion/styled";
import { fontHeaderXxl, Tag } from "@czi-sds/components";

export const TOP_PADDING_PX = 32;
export const LEFT_RIGHT_PADDING_PX = 40;
export const TISSUE_CARD_MAX_WIDTH = 1440;

export const Wrapper = styled.div`
  display: flex;
  flex-direction: column;
  align-self: stretch;
  min-width: 640px;
  max-width: ${TISSUE_CARD_MAX_WIDTH}px;
  margin: auto;
  padding: ${TOP_PADDING_PX}px ${LEFT_RIGHT_PADDING_PX}px 0px
    ${LEFT_RIGHT_PADDING_PX}px;
`;

export const TissueCardHeaderInnerWrapper = styled.div`
  display: flex;
  column-gap: 8px;
  align-items: center;
`;

export const DescriptionWrapper = styled.div`
  margin-bottom: 16px;
`;

export const TissueCardHeader = styled.div`
  display: flex;
  column-gap: 8px;
  flex-direction: row;
  align-items: center;
  justify-content: space-between;
`;

export const TissueCardName = styled.div`
  ${fontHeaderXxl}
  font-weight: 700;
`;

export const StyledTag = styled(Tag)`
  height: 24px;
  margin: 0;
`;

export const SearchBarWrapper = styled.div`
  margin-bottom: 16px;
  margin-top: 16px;
  width: 240px;
`;
