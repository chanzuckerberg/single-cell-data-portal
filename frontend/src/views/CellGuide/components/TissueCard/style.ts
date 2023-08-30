import styled from "@emotion/styled";
import {
  CommonThemeProps,
  fontHeaderXxl,
  getSpaces,
  Tag,
} from "@czi-sds/components";

export const TOP_PADDING_PX = 32;
export const LEFT_RIGHT_PADDING_PX = 40;
export const TISSUE_CARD_MAX_WIDTH = 1440;

interface TissueViewProps extends CommonThemeProps {
  skinnyMode: boolean;
}

export const Wrapper = styled.div`
  display: flex;
  flex-direction: column;
  align-self: stretch;
  margin: auto;
  padding: ${TOP_PADDING_PX}px ${LEFT_RIGHT_PADDING_PX}px 0px
    ${LEFT_RIGHT_PADDING_PX}px;
  margin-bottom: 80px;

  ${(props: TissueViewProps) => {
    const { skinnyMode } = props;
    const spaces = getSpaces(props);
    const space = skinnyMode ? spaces?.l : spaces?.xxl;
    const maxWidth = skinnyMode ? "100vw" : "1440px";

    return `
      max-width: ${maxWidth};
      padding: ${TOP_PADDING_PX}px ${space}px 0px ${space}px;
    `;
  }}
`;

export const TissueCardHeaderInnerWrapper = styled.div`
  display: flex;
  column-gap: 8px;
  align-items: center;
`;

export const DescriptionWrapper = styled.div`
  margin-bottom: 0;
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
