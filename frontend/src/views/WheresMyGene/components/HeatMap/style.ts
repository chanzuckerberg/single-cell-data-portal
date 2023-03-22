import styled from "@emotion/styled";
import { X_AXIS_CHART_HEIGHT_PX, Y_AXIS_CHART_WIDTH_PX } from "./utils";
import { LIGHT_GRAY } from "src/components/common/theme";
import { GENE_INFO_BUTTON_PADDING } from "./components/XAxisChart/style";

export const CHART_LEFT_PADDING = 10;

export const SELECTED_STYLE = {
  backgroundColor: LIGHT_GRAY.D,
  fontWeight: "bold" as never,
  fontFamily: "sans-serif",
  fontSize: 12,
  padding: 4,
};

const Y_AXIS_TOP_PADDING = GENE_INFO_BUTTON_PADDING + 5;

export const Container = styled.div`
  height: 75vh;
  width: 100%;
  overflow: scroll;
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
  width: ${Y_AXIS_CHART_WIDTH_PX + CHART_LEFT_PADDING}px;
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
  padding-left: ${CHART_LEFT_PADDING}px;
  padding-top: 5px;
  left: ${Y_AXIS_CHART_WIDTH_PX}px;
  top: ${X_AXIS_CHART_HEIGHT_PX}px;
`;
