import styled from "@emotion/styled";
import { X_AXIS_CHART_HEIGHT_PX, Y_AXIS_CHART_WIDTH_PX } from "../../utils";

const ECHART_AXIS_LABEL_COLOR_HEX = "#6e7079";
const ECHART_AXIS_LABEL_FONT_SIZE = 12;
const X_AXIS_TITLE_HEIGHT_PX = 40;
const TISSUE_BORDER_WIDTH_PX = 4;

export const X_AXIS_CHART_TOTAL_HEIGHT_PX =
  X_AXIS_CHART_HEIGHT_PX + X_AXIS_TITLE_HEIGHT_PX + TISSUE_BORDER_WIDTH_PX;

export const XAxisContainer = styled.div`
  ${xAxisWidth}
  background-color: white;
  height: ${X_AXIS_CHART_HEIGHT_PX + TISSUE_BORDER_WIDTH_PX}px;
  top: ${X_AXIS_TITLE_HEIGHT_PX + TISSUE_BORDER_WIDTH_PX}px;
  position: absolute;
`;

export const XAxisWrapper = styled.div`
  ${xAxisWidthAndOffset}
  background-color: white;
  height: ${X_AXIS_CHART_TOTAL_HEIGHT_PX}px;
  position: absolute;
  z-index: 2;
`;

export const GeneGroupWrapper = styled.div`
  ${xAxisWidth}

  background-color: white;
  border-bottom: ${TISSUE_BORDER_WIDTH_PX}px solid black;
  height: ${X_AXIS_TITLE_HEIGHT_PX}px;
`;

export const GeneGroupName = styled.div`
  font-size: ${ECHART_AXIS_LABEL_FONT_SIZE}px;
  font-weight: bold;
  top: ${X_AXIS_TITLE_HEIGHT_PX - ECHART_AXIS_LABEL_FONT_SIZE - 6}px;
  position: absolute;
`;

export const MarkerGeneHeader = styled.div`
  font-size: ${ECHART_AXIS_LABEL_FONT_SIZE}px;
  font-weight: normal;
  top: ${X_AXIS_TITLE_HEIGHT_PX - 2 * ECHART_AXIS_LABEL_FONT_SIZE - 10}px;
  position: absolute;
  color: ${ECHART_AXIS_LABEL_COLOR_HEX};
`;

export const MarkerGeneHeaderButton = styled.button`
  font-size: ${ECHART_AXIS_LABEL_FONT_SIZE}px;
  border: none;
  font-weight: normal;
  background-color: white;
  cursor: pointer;
  top: ${X_AXIS_TITLE_HEIGHT_PX - 2 * ECHART_AXIS_LABEL_FONT_SIZE - 10}px;
  right: 0;
  position: absolute;
  color: ${ECHART_AXIS_LABEL_COLOR_HEX};
`;

export const FlexColumn = styled.div`
  display: flex;
  flex-direction: column;
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
  top: ${X_AXIS_TITLE_HEIGHT_PX}px;
  text-align: right;
  left: ${Y_AXIS_CHART_WIDTH_PX - 20}px;
  z-index: 2;
`;

// shift left by the amount the mask is widened due to CellCountLabel
export const MaskWrapper = styled.div`
  display: flex;
  margin-left: -${ECHART_AXIS_LABEL_FONT_SIZE}px;
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
