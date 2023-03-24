import styled from "@emotion/styled";
import { FILTERS_PANEL_EXPANDED_WIDTH_PX } from "src/components/common/SideBar";
import { HEADER_HEIGHT_PX } from "src/components/Header/style";
import {
  CONTENT_WRAPPER_LEFT_RIGHT_PADDING,
  CONTENT_WRAPPER_TOP_BOTTOM_PADDING,
} from "src/components/Layout/style";
import { CELL_INFO_SIDEBAR_WIDTH_PX } from "../CellInfoSideBar/style";
import {
  X_AXIS_CHART_HEIGHT_PX,
  Y_AXIS_CHART_WIDTH_PX,
} from "../HeatMap/utils";
import {
  LEGEND_HEIGHT,
  LEGEND_MARGIN_BOTTOM,
} from "../InfoPanel/components/Legend/style";

export const Header = styled.h1`
  margin-bottom: 12px;
  font-size: 36px;
  font-weight: bold;
`;

const GAP = 10;

const OFFSET_FROM_TOP =
  HEADER_HEIGHT_PX +
  LEGEND_HEIGHT +
  LEGEND_MARGIN_BOTTOM +
  CONTENT_WRAPPER_TOP_BOTTOM_PADDING;

interface WrapperProps {
  isHidden: boolean;
}
// (seve): grid should handle the width for us if we set to 100%, but we
// will need to change upstream styling due to position: absolute
export const Wrapper = styled.div`
  position: absolute;
  margin: 0;
  top: ${OFFSET_FROM_TOP}px;

  display: ${({ isHidden }: WrapperProps) => (isHidden ? "none" : "flex")};

  /* This is to make sure getting started box doesn't overlap FMG panel */
  width: calc(
    100vw -
      ${FILTERS_PANEL_EXPANDED_WIDTH_PX +
      CELL_INFO_SIDEBAR_WIDTH_PX +
      CONTENT_WRAPPER_LEFT_RIGHT_PADDING * 2}px
  );

  height: calc(
    100vh - ${OFFSET_FROM_TOP + CONTENT_WRAPPER_TOP_BOTTOM_PADDING}px
  );
`;

export const Content = styled.div`
  display: flex;
  justify-content: space-between;
`;

export const ColumnOne = styled.div`
  ${isHidden}

  flex-shrink: 0;

  width: ${Y_AXIS_CHART_WIDTH_PX}px;

  z-index: 10;

  display: flex;
  flex-direction: column;
`;

export const ColumnTwo = styled.div`
  flex: 1 0;

  z-index: 10;

  display: flex;
  flex-direction: column;
  justify-content: space-evenly;
`;

export const StyledStepOne = styled.div`
  flex: 1 0;
`;

export const StyledStepTwo = styled.div`
  ${isHidden}
  height: ${X_AXIS_CHART_HEIGHT_PX}px;
  margin-left: ${GAP}px;
`;

export const StyledStepThree = styled.div`
  ${isHidden}
  flex: 1;
  margin: ${GAP}px 0 0 ${GAP}px;
`;

interface IsHiddenProps {
  isHidden?: boolean;
}

function isHidden({ isHidden }: IsHiddenProps) {
  if (!isHidden) return null;

  return "visibility: hidden;";
}
