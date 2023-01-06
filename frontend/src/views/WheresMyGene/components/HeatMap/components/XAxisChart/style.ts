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
  color: ${ECHART_AXIS_LABEL_COLOR_HEX};
  width: ${HEAT_MAP_BASE_CELL_WIDTH_PX}px;
  height: ${X_AXIS_CHART_HEIGHT_PX}px;

  background-color: ${({ active }: { active: boolean }) =>
    active ? SELECTED_STYLE.backgroundColor : "white"};
  font: normal
    ${({ active }: { active: boolean }) =>
      active ? SELECTED_STYLE.fontWeight : "normal"}
    ${SELECTED_STYLE.fontSize}px ${SELECTED_STYLE.fontFamily};
  padding: ${({ active }: { active: boolean }) =>
    active ? SELECTED_STYLE.padding : "unset"};

  text-orientation: sideways;
  writing-mode: vertical-rl;
  text-align: right;
  z-index: 2;
  cursor: pointer;
  user-select: none;
  padding-right: 3px;
  text-overflow: ellipsis;
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
