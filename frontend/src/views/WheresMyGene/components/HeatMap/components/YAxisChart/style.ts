import styled from "@emotion/styled";
import { X_AXIS_CHART_HEIGHT_PX, Y_AXIS_CHART_WIDTH_PX } from "../../utils";

const Y_AXIS_TISSUE_WIDTH_PX = 30;

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
  padding-left: 30x;
`;

export const ResetImageWrapper = styled.div`
  margin-top: 3px;
  cursor: pointer;
`;

function yAxisHeight({ height }: { height: number }) {
  return `
    height: ${height - X_AXIS_CHART_HEIGHT_PX}px;
  `;
}
