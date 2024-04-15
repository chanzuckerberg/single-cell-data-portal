import styled from "@emotion/styled";
import {
  CommonThemeProps,
  fontBodyXs,
  fontHeaderXl,
  fontHeaderXxl,
  Tag,
} from "@czi-sds/components";

import RightSideBar from "src/components/common/RightSideBar";
import { HEADER_HEIGHT_PX } from "src/components/Header/style";
import {
  fontWeightSemibold,
  primary400,
  spacesL,
  spacesM,
  spacesS,
  spacesXxl,
} from "src/common/theme";
import { StyledDiv } from "src/views/WheresMyGeneV2/components/ScreenTint/style";
import { keyframes } from "@emotion/react";
import { DEFAULT_ONTOLOGY_WIDTH } from "../common/OntologyDagView/common/constants";
import { SKINNY_MODE_BREAKPOINT_WIDTH } from "./constants";
import { TOP_PADDING_PX } from "./components/CellGuideCardSidebar/style";
import { Autocomplete } from "@mui/material";

export const SIDEBAR_COLUMN_GAP_PX = 120;
export const CELLGUIDE_CARD_MAX_WIDTH = 1440;
// spacing.xxl and spacing.xl
export const LEFT_RIGHT_PADDING_PX_XXL = 40;

const CELLGUIDE_HEADER_PADDING_BOTTOM_PX = 24;
const CELLGUIDE_HEADER_LINE_HEIGHT_PX = 36;

interface CellGuideViewProps extends CommonThemeProps {
  skinnyMode: boolean;
  topMargin?: number;
}

const AUTOCOMPLETE_HEIGHT = 24;
export const StyledAutocomplete = styled(Autocomplete)`
  height: ${AUTOCOMPLETE_HEIGHT}px;
  width: 135px;
  z-index: 1;
  background-color: white;

  & .MuiAutocomplete-inputRoot {
    padding: 0px;
    height: ${AUTOCOMPLETE_HEIGHT}px;
    display: flex;
    align-items: center;
  }

  & .MuiAutocomplete-inputRoot[class*="MuiOutlinedInput-root"] {
    padding: 0px;
  }

  & .MuiAutocomplete-input {
    padding: 0px;
    height: ${AUTOCOMPLETE_HEIGHT}px;
    line-height: ${AUTOCOMPLETE_HEIGHT - 4}px;
  }

  & .MuiInputLabel-root {
    margin-top: -12px;
    z-index: 0;
  }

  & .MuiInputLabel-root.MuiInputLabel-shrink {
    margin-top: 0px;
  }
  .MuiOutlinedInput-input.MuiOutlinedInput-input {
    padding: 0px 8px;
  }
  .MuiAutocomplete-popupIndicator,
  .MuiAutocomplete-clearIndicator {
    width: ${AUTOCOMPLETE_HEIGHT - 4}px;
    height: ${AUTOCOMPLETE_HEIGHT - 4}px;
    margin-top: 4px;
  }
`;

export const CellGuideView = styled.div<CellGuideViewProps>`
  display: flex;
  flex-direction: row;
  column-gap: ${SIDEBAR_COLUMN_GAP_PX}px;
  max-width: 100vw;
  height: auto;

  ${(props) => {
    const { skinnyMode } = props;
    const space = skinnyMode ? spacesL(props) : spacesXxl(props);
    return `
      padding: 0px ${space}px;
    `;
  }}
`;

export const CellGuideWrapper = styled.div<CellGuideViewProps>`
  margin: 0 auto 80px;
  max-width: ${(props) =>
    props.skinnyMode
      ? `${DEFAULT_ONTOLOGY_WIDTH}px`
      : `${CELLGUIDE_CARD_MAX_WIDTH}px`};
`;

export const Wrapper = styled.div<CellGuideViewProps>`
  display: flex;
  flex-direction: column;
  align-self: stretch;
  overflow-x: hidden;

  ${(props) => {
    const { skinnyMode, topMargin } = props;
    const maxWidth = skinnyMode
      ? `${DEFAULT_ONTOLOGY_WIDTH}px`
      : `${SKINNY_MODE_BREAKPOINT_WIDTH + SIDEBAR_COLUMN_GAP_PX}px`;

    const marginTop =
      topMargin && !skinnyMode
        ? `${topMargin}px`
        : `${
            CELLGUIDE_HEADER_LINE_HEIGHT_PX +
            CELLGUIDE_HEADER_PADDING_BOTTOM_PX +
            TOP_PADDING_PX
          }px`;
    return `
    max-width: ${maxWidth};
    margin-top: ${marginTop};
    `;
  }}
`;

export const NavBarDropdownWrapper = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: flex-end;
  padding-bottom: ${spacesM}px;
`;

export const CellGuideCardHeaderInnerWrapper = styled.div`
  display: flex;
  column-gap: ${spacesS}px;
  align-items: center;
`;

interface CellGuideCardHeaderProps extends CommonThemeProps {
  width: number;
}

export const CellGuideCardHeader = styled.div<CellGuideCardHeaderProps>`
  display: flex;
  column-gap: ${spacesS}px;
  flex-direction: row;
  align-items: center;
  justify-content: space-between;
  position: fixed;
  top: ${HEADER_HEIGHT_PX}px;
  padding-bottom: ${CELLGUIDE_HEADER_PADDING_BOTTOM_PX}px;
  line-height: ${CELLGUIDE_HEADER_LINE_HEIGHT_PX}px;
  padding-top: ${TOP_PADDING_PX}px;
  padding-left: ${spacesXxl}px;
  background: linear-gradient(180deg, #fff 89.82%, rgba(255, 255, 255, 0) 100%);
  z-index: 10;
  width: ${(props) => `${props.width + LEFT_RIGHT_PADDING_PX_XXL}px`};
  @media (min-width: ${CELLGUIDE_CARD_MAX_WIDTH}px) {
    left: calc(50vw - ${CELLGUIDE_CARD_MAX_WIDTH / 2}px);
  }
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

export const StyledCellTagSideBar = styled(Tag)`
  .MuiChip-label {
    color: ${primary400};
    font-weight: 500 !important;
    ${fontBodyXs}
  }
  background-color: #e0f0ff;
  &:hover {
    cursor: default;
    background-color: #e0f0ff;
    .MuiChip-label {
      color: ${primary400};
    }
  }
`;

export const StyledGeneTagSideBar = styled(Tag)`
  .MuiChip-label {
    color: #8f5aff;
    font-weight: 500 !important;
    ${fontBodyXs}
  }
  background-color: #8f5aff26;
  &:hover {
    cursor: default;
    background-color: #8f5aff26;
    .MuiChip-label {
      color: #8f5aff;
    }
  }
`;

export const SearchBarPositioner = styled.div`
  display: flex;
  justify-content: flex-end;
`;

export const SearchBarWrapper = styled.div`
  margin-bottom: 20px;
  width: 100%;
`;

const slideIn = keyframes`
  from {
    transform: translateX(100%);
  }
  to {
    transform: translateX(0);
  }
`;

interface StyledRightSideBarProps extends CommonThemeProps {
  skinnyMode: boolean;
}

export const StyledRightSideBar = styled(RightSideBar)<StyledRightSideBarProps>`
  position: fixed;
  right: 0;
  height: 100vh;
  background-color: white;
  z-index: 10;
  top: ${HEADER_HEIGHT_PX}px;
  animation: ${slideIn} 0.2s ease-in-out forwards;
  ${(props) => {
    if (props.skinnyMode) {
      return `
        width: 100vw;
        height: calc(100vh - ${HEADER_HEIGHT_PX}px);
      `;
    }
  }}
`;

export const MobileSearchTint = styled(StyledDiv)`
  z-index: 1;
  background: rgba(0, 0, 0, 0.3);
  position: fixed;
`;

export const MobileTooltipTitle = styled.div`
  ${fontHeaderXl}
  font-weight: ${fontWeightSemibold};
`;

export const MobileTooltipWrapper = styled.div`
  top: ${HEADER_HEIGHT_PX}px;
  position: fixed;
  z-index: 99;
  background-color: white;
  display: flex;
  flex-direction: column;
  height: calc(100vh - ${HEADER_HEIGHT_PX}px);
  overscroll-behavior: none;
  width: 100vw;
  padding: ${spacesL}px;
  gap: ${spacesL}px;
  overflow: auto;
`;

export const MobileTooltipHeader = styled.div`
  display: flex;
  justify-content: space-between;
  align-items: center;
`;
