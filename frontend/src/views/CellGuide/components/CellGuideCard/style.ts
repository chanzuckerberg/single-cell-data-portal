import styled from "@emotion/styled";
import {
  CommonThemeProps,
  fontHeaderXl,
  fontHeaderXxl,
  getSpaces,
  Tag,
} from "@czi-sds/components";

import RightSideBar from "src/components/common/RightSideBar";
import { HEADER_HEIGHT_PX } from "src/components/Header/style";
import Synonyms from "src/components/Synonyms";
import { spacesL, spacesXxs } from "src/common/theme";
import { StyledDiv } from "src/views/WheresMyGene/components/ScreenTint/style";
import OntologyId from "src/components/OntologyId";

export const TOP_PADDING_PX = 32;
export const SIDEBAR_COLUMN_GAP_PX = 120;

// spacing.xxl and spacing.xl
export const LEFT_RIGHT_PADDING_PX = 40;
export const LEFT_RIGHT_PADDING_PX_SKINNY_MODE = 24;

interface CellGuideViewProps extends CommonThemeProps {
  skinnyMode: boolean;
}

export const CellGuideView = styled.div`
  display: flex;
  flex-direction: row;
  column-gap: ${SIDEBAR_COLUMN_GAP_PX}px;
  margin: auto;
  margin-bottom: 80px;

  ${(props: CellGuideViewProps) => {
    const { skinnyMode } = props;
    const spaces = getSpaces(props);
    const space = skinnyMode ? spaces?.l : spaces?.xxl;
    const maxWidth = skinnyMode ? "100vw" : "1440px";

    return `
    max-width: ${maxWidth};
    padding: ${TOP_PADDING_PX}px ${space}px 0px
      ${space}px;

    `;
  }}
`;

export const Wrapper = styled.div`
  display: flex;
  flex-direction: column;
  align-self: stretch;
  margin-bottom: 80px;
`;

export const CellGuideCardHeaderInnerWrapper = styled.div`
  display: flex;
  column-gap: 8px;
  align-items: center;
`;

export const CellGuideCardHeader = styled.div`
  display: flex;
  column-gap: 8px;
  flex-direction: row;
  align-items: center;
  justify-content: space-between;
`;

export const CellGuideCardName = styled.h1`
  ${fontHeaderXxl}
  font-weight: 700;
  margin-bottom: 0;
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
  width: 100%;
`;

export const StyledRightSideBar = styled(RightSideBar)`
  position: fixed;
  right: 0;
  height: 100vh;
  background-color: white;
  z-index: 10; // Must be over mobile cell guide nav bar
  top: ${HEADER_HEIGHT_PX}px;

  ${(props: { skinnyMode: boolean }) => {
    if (props.skinnyMode) {
      return `
        width: 100vw;
        height: calc(100vh - ${HEADER_HEIGHT_PX}px);
      `;
    }
  }}
`;

export const StyledSynonyms = styled(Synonyms)`
  margin-top: ${spacesXxs}px;
  margin-left: ${spacesL}px;
`;

export const StyledOntologyId = styled(OntologyId)`
  margin-top: ${spacesXxs}px;
  margin-left: ${spacesL}px;
`;

export const MobileSearchTint = styled(StyledDiv)`
  z-index: 1;
  background: rgba(0, 0, 0, 0.3);
  position: fixed;
`;

export const FlexContainer = styled.div`
  display: flex;
  flex-direction: column;
`;

export const MobileTooltipTitle = styled.div`
  ${fontHeaderXl}
  font-weight: 600;
`;

export const MobileTooltipWrapper = styled.div`
  top: ${HEADER_HEIGHT_PX}px;
  position: fixed;
  z-index: 99;
  background-color: white;
  display: flex;
  flex-direction: column;
  height: calc(100vh - ${HEADER_HEIGHT_PX}px);
  width: 100vw;
  padding: ${spacesL}px;
  gap: ${spacesL}px;
`;

export const MobileTooltipHeader = styled.div`
  display: flex;
  justify-content: space-between;
  align-items: center;
`;
