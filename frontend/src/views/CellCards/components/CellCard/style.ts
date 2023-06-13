import styled from "@emotion/styled";
import {
  Button,
  CommonThemeProps,
  fontBodyS,
  fontHeaderXxl,
  getSpaces,
  Tag,
} from "@czi-sds/components";

import RightSideBar from "src/components/common/RightSideBar";
import { HEADER_HEIGHT_PX } from "src/components/Header/style";

export const TOP_PADDING_PX = 32;
export const SIDEBAR_COLUMN_GAP_PX = 120;

// spacing.xxl and spacing.xl
export const LEFT_RIGHT_PADDING_PX = 40;
export const LEFT_RIGHT_PADDING_PX_SKINNY_MODE = 24;

interface CellCardsViewProps extends CommonThemeProps {
  skinnyMode: boolean;
}

export const CellCardsView = styled.div`
  display: flex;
  flex-direction: row;
  column-gap: ${SIDEBAR_COLUMN_GAP_PX}px;
  margin: auto;
  max-width: 1440px;

  ${(props: CellCardsViewProps) => {
    const { skinnyMode } = props;
    const spaces = getSpaces(props);
    const space = skinnyMode ? spaces?.xl : spaces?.xxl;

    return `
    padding: ${TOP_PADDING_PX}px ${space}px 0px
      ${space}px;

    `;
  }}
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

export const SearchBarPositioner = styled.div`
  display: flex;
  justify-content: flex-end;
`;

export const SearchBarWrapper = styled.div`
  margin-bottom: 20px;
  width: 240px;
`;

export const SuggestChangeButton = styled(Button)`
  ${fontBodyS}
  font-weight: 500;
  text-transform: capitalize;
`;

export const StyledRightSideBar = styled(RightSideBar)`
  position: fixed;
  top: ${HEADER_HEIGHT_PX}px;
  right: 0;
  height: 100vh;
  background-color: white;
`;
