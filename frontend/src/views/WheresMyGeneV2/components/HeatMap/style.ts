import styled from "@emotion/styled";
import { CommonThemeProps, TagFilter, fontHeaderXl } from "@czi-sds/components";

import { gray100, gray400, gray300, spacesM, spacesS } from "src/common/theme";
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
  visible: boolean;
}

export const ChartWrapper = styled.div<ChartWrapperProps>`
  position: absolute;
  padding-left: ${CHART_PADDING_PX + GENE_CHART_LEFT_OFFSET_PX}px;
  padding-right: ${CHART_PADDING_PX}px;
  padding-top: ${(props) => xAxisOffset(props) + PADDING_UNDER_HEADERS_PX}px;
  left: ${Y_AXIS_CHART_WIDTH_PX}px;
  top: ${(props) => props.top}px;
  visibility: ${(props) => (props.visible ? "visible" : "hidden")};
  width: 100%;
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

interface LoadingContainerProps extends CommonThemeProps {
  height: number;
  width: number;
}

export const LoadingContainer = styled.div<LoadingContainerProps>`
  position: absolute;
  display: flex;
  width: ${(props) => props.width * 20}px;
  height: ${(props) => props.height * 20}px;
  background-color: ${gray100};
`;

interface LoadingProps extends CommonThemeProps {
  left: number;
}

export const LoadingWrapper = styled.div<LoadingProps>`
  display: flex;
  flex-direction: column;
  justify-content: center;
  align-items: center;
  position: fixed;
  top: 600px;
  left: ${(props) =>
    (props.left < 13 ? 640 : 510) + (props.left * 19 * 1.27) / 2}px;
`;

interface LoadingLabelProps extends CommonThemeProps {
  hidden: boolean;
}

export const LoadingLabel = styled.div<LoadingLabelProps>`
  color: black;
  ${fontHeaderXl}
  padding-top: 20px;
  visibility: ${(props) => (props.hidden ? "hidden" : "visible")};
`;

export const LoadingSpinner = styled.div`
  --d: 10.6px;
  width: 1.9px;
  height: 1.9px;
  border-radius: 50%;
  color: ${gray400};
  box-shadow:
    calc(1 * var(--d)) calc(0 * var(--d)) 0 0,
    calc(0.707 * var(--d)) calc(0.707 * var(--d)) 0 0.5px,
    calc(0 * var(--d)) calc(1 * var(--d)) 0 1px,
    calc(-0.707 * var(--d)) calc(0.707 * var(--d)) 0 1.4px,
    calc(-1 * var(--d)) calc(0 * var(--d)) 0 1.9px,
    calc(-0.707 * var(--d)) calc(-0.707 * var(--d)) 0 2.4px;
  animation: spinner 1s infinite steps(8);

  @keyframes spinner {
    100% {
      transform: rotate(1turn);
    }
  }
`;
