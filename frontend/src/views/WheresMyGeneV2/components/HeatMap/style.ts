import styled from "@emotion/styled";
import { gray300 } from "src/common/theme";
import { X_AXIS_CHART_HEIGHT_PX } from "src/views/WheresMyGeneV2/common/constants";
import {
  CELL_TYPE_FILTER_WIDTH_PX,
  DIVIDER_MARGIN_PX,
  Y_AXIS_CHART_WIDTH_PX,
} from "src/views/WheresMyGene/components/HeatMap/utils";

enum ZIndex {
  XAxisWrapper = 2,
  TopLeftCornerMask = 3,
  Divider = 4,
}

export const CellTypeFilterContainer = styled.div`
  height: 100%;
  width: ${CELL_TYPE_FILTER_WIDTH_PX}px;
`;

export const Divider = styled.div`
  height: 100%;
  position: absolute;
  left: ${CELL_TYPE_FILTER_WIDTH_PX + DIVIDER_MARGIN_PX}px;
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

interface TopLeftCornerMaskProps {
  height: number;
}

export const TopLeftCornerMask = styled.div<TopLeftCornerMaskProps>`
  position: absolute;
  background-color: white;
  z-index: ${ZIndex.TopLeftCornerMask};
  top: 0px;
  left: 0px;
  width: ${Y_AXIS_CHART_WIDTH_PX}px;
  height: ${(props) => props.height || X_AXIS_CHART_HEIGHT_PX}px;
  min-height: ${(props) => props.height}px;
  display: flex;
  flex-direction: row;
  justify-content: space-between;
  align-items: end;
`;
