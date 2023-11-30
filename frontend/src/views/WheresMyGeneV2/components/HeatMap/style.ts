import styled from "@emotion/styled";
import { CommonThemeProps, TagFilter } from "@czi-sds/components";

import { gray300, spacesM, spacesS } from "src/common/theme";
import { HEADER_HEIGHT_PX } from "src/components/Header/style";
import {
  CONTENT_WRAPPER_LEFT_RIGHT_PADDING_PX,
  CONTENT_WRAPPER_TOP_BOTTOM_PADDING_PX,
} from "src/components/Layout/style";
import { X_AXIS_CHART_HEIGHT_PX } from "src/views/WheresMyGeneV2/common/constants";
import {
  CELL_TYPE_FILTER_WIDTH_PX,
  DIVIDER_LEFT_POSITION_PX,
  GENE_CHART_LEFT_OFFSET_PX,
  Y_AXIS_CHART_WIDTH_PX,
} from "src/views/WheresMyGeneV2/components/HeatMap/utils";
import { LEGEND_HEIGHT_PX } from "../InfoPanel/components/Legend/style";
import { LEGEND_MARGIN_BOTTOM_PX } from "../../style";
import { LIGHT_GRAY } from "src/components/common/theme";
import { Skeleton } from "@mui/material";

export function xAxisOffset(props: CommonThemeProps) {
  /**
   * (thuang): This offset is to make sure the x-axis label doesn't overlap the
   * gene search bar.
   */
  return spacesM?.(props) || 0;
}

const PADDING_UNDER_HEADERS_PX = 5;
const TOP_LEFT_CORNER_OFFSET_PX = 6;

export const SELECTED_STYLE = {
  backgroundColor: LIGHT_GRAY.D,
  fontWeight: "bold" as never,
  fontFamily: "sans-serif",
  fontSize: 12,
  padding: 4,
};

enum ZIndex {
  XAxisWrapper = 2,
  TopLeftCornerMask = 3,
  Divider = 4,
}

export const CHART_PADDING_PX = 10;

export const CellTypeFilterContainer = styled.div`
  height: 100%;
  width: ${CELL_TYPE_FILTER_WIDTH_PX}px;
`;

const CELL_TYPE_SEARCH_BOX_HEIGHT_PX = 37;

export const CellTypeTagContainer = styled.div`
  overflow-y: auto;
  height: calc(100% - ${CELL_TYPE_SEARCH_BOX_HEIGHT_PX}px);
  padding: ${spacesS}px;
`;

export const Container = styled.div`
  /* Fixes a bug where PNG was getting cut off */
  width: calc(100% + ${CONTENT_WRAPPER_LEFT_RIGHT_PADDING_PX}px);

  height: calc(
    100vh -
      ${HEADER_HEIGHT_PX +
      LEGEND_HEIGHT_PX +
      LEGEND_MARGIN_BOTTOM_PX +
      CONTENT_WRAPPER_TOP_BOTTOM_PADDING_PX}px
  );
  overflow: auto;
  position: relative;
`;

export const ContainerWrapper = styled.div`
  position: relative;
`;

export const Divider = styled.div`
  height: 100%;
  position: absolute;
  left: ${DIVIDER_LEFT_POSITION_PX}px;
  width: 1px;
  top: 0;
  border-right: solid 0.5px ${gray300};
  z-index: ${ZIndex.Divider};
`;

export const XAxisWrapper = styled.div`
  display: flex;
  background-color: white;
  flex-direction: row;
  position: sticky;
  top: 0;
  width: 100%;
  z-index: ${ZIndex.XAxisWrapper};
`;

interface TopLeftCornerMaskProps extends CommonThemeProps {
  height: number;
}

export const TopLeftCornerMask = styled.div<TopLeftCornerMaskProps>`
  position: absolute;
  background-color: white;
  z-index: ${ZIndex.TopLeftCornerMask};
  top: 0px;
  left: 0px;
  width: ${Y_AXIS_CHART_WIDTH_PX + TOP_LEFT_CORNER_OFFSET_PX}px;
  height: ${(props) =>
    props.height + xAxisOffset(props) || X_AXIS_CHART_HEIGHT_PX}px;
  min-height: ${(props) => props.height}px;
  display: flex;
  flex-direction: row;
  justify-content: space-between;
  align-items: end;
`;

interface ChartWrapperProps extends CommonThemeProps {
  top: number;
  hidden: boolean;
}

export const ChartWrapper = styled.div<ChartWrapperProps>`
  position: absolute;
  padding-left: ${CHART_PADDING_PX + GENE_CHART_LEFT_OFFSET_PX}px;
  padding-right: ${CHART_PADDING_PX}px;
  padding-top: ${(props) => xAxisOffset(props) + PADDING_UNDER_HEADERS_PX}px;
  left: ${Y_AXIS_CHART_WIDTH_PX}px;
  top: ${(props) => props.top}px;
  visibility: ${(props) => (props.hidden ? "hidden" : "visible")};
`;

export const StyledTag = styled(TagFilter)`
  max-width: 258px;
`;

interface XAxisMaskProps {
  height: number;
}

export const XAxisMask = styled.div<XAxisMaskProps>`
  width: ${Y_AXIS_CHART_WIDTH_PX + CHART_PADDING_PX}px;
  height: ${(props) => props.height}px;
`;

interface YAxisWrapperProps extends CommonThemeProps {
  top: number;
}

export const YAxisWrapper = styled.div<YAxisWrapperProps>`
  width: ${Y_AXIS_CHART_WIDTH_PX}px;
  position: sticky;
  top: ${(props) => props.top ?? X_AXIS_CHART_HEIGHT_PX}px;
  left: 0;
  z-index: 1;
  padding-top: ${(props) => xAxisOffset(props) + PADDING_UNDER_HEADERS_PX}px;
  /* Somehow Firefox requires this to scroll */
  overflow: hidden;
`;

export const StyledSkeleton = styled(Skeleton)`
  margin: 0px 0px 4px 1px;
`;

export const SkeletonContainer = styled.div`
  display: flex;
  flex-direction: column;
`;

export const SkeletonWrapper = styled.div`
  display: flex;
  flex-direction: row;
`;
