import styled from "@emotion/styled";
import { SELECTED_STYLE } from "../../style";
import {
  HEAT_MAP_BASE_CELL_WIDTH_PX,
  Y_AXIS_CHART_WIDTH_PX,
} from "../../utils";

export const ECHART_AXIS_LABEL_COLOR_HEX = "#6e7079";
export const ECHART_AXIS_LABEL_FONT_SIZE_PX = 12;
export const GENE_INFO_BUTTON_PADDING_PX = 12;

interface XAxisContainerProps {
  height: number;
  width: number;
}

export const XAxisContainer = styled.div<XAxisContainerProps>`
  ${xAxisWidth}
  background-color: white;
  height: ${(props) => props.height}px;
  top: 0px;
  position: absolute;
  display: flex;
  flex-direction: row;
  justify-content: space-between;
`;

interface XAxisWrapperProps {
  height: number;
  width: number;
  left: number;
}

export const XAxisWrapper = styled.div<XAxisWrapperProps>`
  ${xAxisWidthAndOffset}
  background-color: white;
  height: ${(props) => props.height}px;
  position: absolute;
  z-index: 2;
`;

export const XAxisLabel = styled.div`
  height: 100%;
  width: ${HEAT_MAP_BASE_CELL_WIDTH_PX}px;
  writing-mode: vertical-rl;
  color: ${ECHART_AXIS_LABEL_COLOR_HEX};
  display: flex;
  justify-content: end;
  align-items: center;
  margin-bottom: 2px;
`;

export const XAxisGeneName = styled.span`
  transform: scale(-1, -1);
  ${selectedStyle}
`;

export const InfoButtonWrapper = styled.div`
  transform: scale(1, 1);
  cursor: pointer;
  display: flex;
  margin-bottom: 4px;
  margin-top: 4px;
  justify-content: center;
`;

export const DeleteButtonWrapper = styled.div`
  cursor: pointer;
`;

export const GeneButtonStyle = styled.div`
  background-color: white;
  border: none;
  z-index: 2;
  display: inline-flex;
  justify-content: space-between;
  white-space: nowrap;
  overflow: hidden;
  flex-direction: column;
  align-items: center;

  :hover {
    #gene-hover-container {
      visibility: visible;
    }
  }
`;

export const HoverContainer = styled.div`
  display: flex;
  flex-direction: column;
  visibility: hidden;
  position: absolute;
`;

export const CellCountLabel = styled.div`
  font: ${ECHART_AXIS_LABEL_FONT_SIZE_PX}px sans-serif;
  color: ${ECHART_AXIS_LABEL_COLOR_HEX};
  width: ${HEAT_MAP_BASE_CELL_WIDTH_PX}px;
  height: 100%;
  writing-mode: vertical-rl;
  position: absolute;
  left: ${Y_AXIS_CHART_WIDTH_PX - HEAT_MAP_BASE_CELL_WIDTH_PX}px;
  transform: scale(-1, -1);
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
