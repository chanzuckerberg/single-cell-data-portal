import styled from "@emotion/styled";
import { FILTERS_PANEL_EXPANDED_WIDTH_PX } from "src/components/common/SideBar";
import { HEADER_HEIGHT_PX } from "src/components/Header/style";
import {
  CONTENT_WRAPPER_LEFT_RIGHT_PADDING_PX,
  CONTENT_WRAPPER_TOP_BOTTOM_PADDING_PX,
} from "src/components/Layout/style";
import { LEGEND_MARGIN_BOTTOM_PX } from "../../style";
import { CELL_INFO_SIDEBAR_WIDTH_PX } from "../CellInfoSideBar/style";
import { Y_AXIS_CHART_WIDTH_PX } from "../HeatMap/utils";
import { LEGEND_HEIGHT_PX } from "../InfoPanel/components/Legend/style";

export const Header = styled.h1`
  margin-bottom: 12px;
  font-size: 36px;
  font-weight: bold;
`;

// Matches padding/gap between heat map chart and x/y axis components
const GAP_PX = 5;

const Z_INDEX = 10;

const OFFSET_FROM_TOP_PX =
  HEADER_HEIGHT_PX +
  LEGEND_HEIGHT_PX +
  LEGEND_MARGIN_BOTTOM_PX +
  CONTENT_WRAPPER_TOP_BOTTOM_PADDING_PX;

interface WrapperProps {
  isHidden: boolean;
}
// (seve): grid should handle the width for us if we set to 100%, but we
// will need to change upstream styling due to position: absolute
export const Wrapper = styled.div`
  position: absolute;
  margin: 0;
  top: ${OFFSET_FROM_TOP_PX}px;

  display: ${({ isHidden }: WrapperProps) => (isHidden ? "none" : "flex")};

  /* This is to make sure getting started box doesn't overlap FMG panel */
  width: calc(
    100vw -
      ${FILTERS_PANEL_EXPANDED_WIDTH_PX +
      CELL_INFO_SIDEBAR_WIDTH_PX +
      CONTENT_WRAPPER_LEFT_RIGHT_PADDING_PX * 2}px
  );

  height: calc(
    100vh - ${OFFSET_FROM_TOP_PX + CONTENT_WRAPPER_TOP_BOTTOM_PADDING_PX}px
  );
`;

export const Content = styled.div`
  display: flex;
  justify-content: space-between;
`;

export const ColumnOne = styled.div`
  flex-shrink: 0;

  width: ${Y_AXIS_CHART_WIDTH_PX}px;

  display: flex;
  flex-direction: column;
`;

export const ColumnTwo = styled.div`
  flex: 1 0;

  display: flex;
  flex-direction: column;
  justify-content: space-evenly;
`;

export const StyledStepOne = styled.div`
  ${isHidden}
  z-index: ${Z_INDEX};
  flex: 1;
`;

type Props = {
  minHeight: number;
  isHidden: boolean;
};

export const StyledStepTwo = styled.div<Props>`
  ${isHidden}
  z-index: ${Z_INDEX};
  flex: 0;
  min-height: ${(props) => props.minHeight}px;
  margin-left: ${GAP_PX}px;
`;

export const StyledStepThree = styled.div`
  ${isHidden}
  z-index: ${Z_INDEX};
  flex: 1;

  margin-left: ${GAP_PX}px;
  margin-top: ${GAP_PX}px;
`;

interface IsHiddenProps {
  isHidden?: boolean;
}

function isHidden({ isHidden }: IsHiddenProps) {
  if (!isHidden) return null;

  return "visibility: hidden;";
}
