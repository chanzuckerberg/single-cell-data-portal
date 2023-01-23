import styled from "@emotion/styled";
import { SELECTED_STYLE } from "../../style";
import { X_AXIS_CHART_HEIGHT_PX, Y_AXIS_CHART_WIDTH_PX } from "../../utils";
import { HEAT_MAP_BASE_CELL_WIDTH_PX } from "../../utils";

export const ECHART_AXIS_LABEL_COLOR_HEX = "#6e7079";
const ECHART_AXIS_LABEL_FONT_SIZE = 12;

export const XAxisContainer = styled.div`
  ${xAxisWidth}
  background-color: white;
  height: ${X_AXIS_CHART_HEIGHT_PX}px;
  top: 0px;
  position: absolute;
  display: flex;
  flex-direction: row;
  justify-content: space-between;
`;

export const XAxisWrapper = styled.div`
  ${xAxisWidthAndOffset}
  background-color: white;
  height: ${X_AXIS_CHART_HEIGHT_PX}px;
  position: absolute;
  z-index: 2;
`;

export const XAxisLabel = styled.div`
  ${selectedStyle}
  width: ${HEAT_MAP_BASE_CELL_WIDTH_PX}px;
  text-orientation: sideways;
  writing-mode: vertical-rl;
  display: inline-block;
  user-select: none;
  color: ${ECHART_AXIS_LABEL_COLOR_HEX};
`;

export const GeneButtonStyle = styled.button`
  cursor: pointer;
  background-color: white;
  border: none;
  z-index: 2;
  display: inline-flex;
  justify-content: center;
  align-items: end;
  white-space: nowrap;
  overflow: hidden;
`;

// adjust the left position of CellCountLabel by -20 to center it properly
export const CellCountLabel = styled.div`
  font: ${ECHART_AXIS_LABEL_FONT_SIZE}px sans-serif;
  color: ${ECHART_AXIS_LABEL_COLOR_HEX};
  width: ${ECHART_AXIS_LABEL_FONT_SIZE}px;
  height: ${X_AXIS_CHART_HEIGHT_PX}px;
  background-color: white;
  text-orientation: sideways;
  writing-mode: vertical-rl;
  padding-top: 16px;
  position: absolute;
  top: 0px;
  text-align: right;
  left: ${Y_AXIS_CHART_WIDTH_PX - 20}px;
  z-index: 2;
`;

function xAxisWidthAndOffset({ width, left }: { width: number; left: number }) {
  return `
    width: ${width}px;
    left: ${left}px;
  `;
}

function xAxisWidth({ width }: { width: number }) {
  return `
    width: ${width}px;
  `;
}

function selectedStyle({ font, active }: { font: string; active: boolean }) {
  return `
    font: ${font};
    background-color: ${active ? SELECTED_STYLE.backgroundColor : "white"};
  `;
}
