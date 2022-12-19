import Image from "next/image";
import styled from "@emotion/styled";
import { SELECTED_STYLE, HEAT_MAP_BASE_CELL_PX, X_AXIS_CHART_HEIGHT_PX, Y_AXIS_CHART_WIDTH_PX } from "../../utils";

export const Y_AXIS_TISSUE_WIDTH_PX = 30;

export const Wrapper = styled.div`
  display: flex;
  margin-bottom: ${X_AXIS_CHART_HEIGHT_PX}px;
  margin-right: ${Y_AXIS_TISSUE_WIDTH_PX}px;
  width: ${Y_AXIS_CHART_WIDTH_PX}px;
`;

const TISSUE_BORDER_WIDTH_PX = 5;

export const TissueWrapper = styled.div`
  ${yAxisHeight}

  background-color: white;
  border-right: ${TISSUE_BORDER_WIDTH_PX}px solid black;
  width: ${Y_AXIS_TISSUE_WIDTH_PX}px;
  padding-left: ${TISSUE_BORDER_WIDTH_PX}px;
`;

export const TissueName = styled.div`
  text-orientation: sideways;
  writing-mode: vertical-rl;
  transform: rotate(180deg);
  font-size: 12px;
  font-weight: bold;
`;

export const Container = styled.div`
  ${yAxisHeight}

  background-color: white;
  width: ${Y_AXIS_CHART_WIDTH_PX}px;
`;

export const CellTypeButtonStyle = styled.button`
  height: ${HEAT_MAP_BASE_CELL_PX}px;
  background-color: ${({active}: {active: boolean})=> active ? SELECTED_STYLE.backgroundColor : "white"};
  font: normal ${({active}: {active: boolean})=> active ? SELECTED_STYLE.fontWeight : "normal"} ${SELECTED_STYLE.fontSize}px ${SELECTED_STYLE.fontFamily};
  padding: ${({active}: {active: boolean})=> active ? SELECTED_STYLE.padding : "unset"};
  white-space: pre;
  cursor: pointer;
  border: none;
  width: 100%;
  color: #6E7079;
  text-align: left;
`

export const CellCountLabelStyle = styled.div`
  height: ${HEAT_MAP_BASE_CELL_PX}px;
  background-color: white;
  font: normal normal ${SELECTED_STYLE.fontSize}px ${SELECTED_STYLE.fontFamily};
  white-space: pre;
  border: none;
  color: #6E7079;
  text-align: right;
  padding-top: 3px;
`
export const FlexRowJustified = styled.div`
  display: flex;
  flex-direction: row;
  justify-content: space-between;
  padding-left: 5px;
  width: 100%;
`

export const FlexRow = styled.div`
  display: flex;
  flex-direction: row;
`

export const InfoButtonWrapper = styled.div`
  padding-left: 2px;
  cursor: pointer;
`;

export const ResetImageWrapper = styled.div`
  margin-top: 3px;
  cursor: pointer;
`;

function yAxisHeight({ height }: { height: number }) {
  return `
    height: ${height}px;
  `;
}

export const StyledImage = styled(Image)`
  :hover {
    filter: brightness(0);
  }
`;
