import styled from "@emotion/styled";
import { X_AXIS_CHART_HEIGHT_PX, Y_AXIS_CHART_WIDTH_PX } from "./utils";
import { LIGHT_GRAY } from "src/components/common/theme";
import { LEGEND_HEIGHT_PX } from "../InfoPanel/components/Legend/style";
import { HEADER_HEIGHT_PX } from "src/components/Header/style";
import {
  CONTENT_WRAPPER_LEFT_RIGHT_PADDING_PX,
  CONTENT_WRAPPER_TOP_BOTTOM_PADDING_PX,
} from "src/components/Layout/style";
import { LEGEND_MARGIN_BOTTOM_PX } from "../../style";

export const CHART_PADDING_PX = 10;

export const SELECTED_STYLE = {
  backgroundColor: LIGHT_GRAY.D,
  fontWeight: "bold" as never,
  fontFamily: "sans-serif",
  fontSize: 12,
  padding: 4,
};

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

export const TopLeftCornerMask = styled.div`
  position: absolute;
  background-color: white;
  z-index: 3;
  top: 0px;
  left: 0px;
  width: ${Y_AXIS_CHART_WIDTH_PX}px;
  height: ${X_AXIS_CHART_HEIGHT_PX}px;
`;

export const YAxisWrapper = styled.div`
  width: ${Y_AXIS_CHART_WIDTH_PX}px;
  position: sticky;
  top: ${X_AXIS_CHART_HEIGHT_PX}px;
  left: 0;
  z-index: 1;
  padding-top: 5px;
  /* Somehow Firefox requires this to scroll */
  overflow: hidden;
`;

export const XAxisMask = styled.div`
  width: ${Y_AXIS_CHART_WIDTH_PX + CHART_PADDING_PX}px;
  height: ${X_AXIS_CHART_HEIGHT_PX}px;
`;

export const XAxisWrapper = styled.div`
  display: flex;
  background-color: white;
  flex-direction: row;
  position: sticky;
  top: 0;
  width: 100%;
  z-index: 2;
`;

export const ChartWrapper = styled.div`
  position: absolute;
  padding-left: ${CHART_PADDING_PX}px;
  padding-right: ${CHART_PADDING_PX}px;
  padding-top: 5px;
  left: ${Y_AXIS_CHART_WIDTH_PX}px;
  top: ${X_AXIS_CHART_HEIGHT_PX}px;
`;
