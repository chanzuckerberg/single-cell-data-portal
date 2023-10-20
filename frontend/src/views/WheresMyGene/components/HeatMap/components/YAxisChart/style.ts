import Image from "next/image";
import styled from "@emotion/styled";
import { HEAT_MAP_BASE_CELL_PX, Y_AXIS_CHART_WIDTH_PX } from "../../utils";
import { ECHART_AXIS_LABEL_COLOR_HEX } from "../XAxisChart/style";
import { SELECTED_STYLE } from "../../style";
import { MARGIN_BETWEEN_HEATMAPS } from "src/views/WheresMyGene/common/constants";
import { fontWeightSemibold, gray500 } from "src/common/theme";
import { CommonThemeProps } from "@czi-sds/components";

export const Y_AXIS_TISSUE_WIDTH_PX = 30;

export const Wrapper = styled.div`
  display: flex;
  margin-bottom: ${MARGIN_BETWEEN_HEATMAPS}px;
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

  /* Fixes bug where some tissue names were wrapping on the last character for png download */
  white-space: nowrap;
`;

export const Container = styled.div`
  ${yAxisHeight}

  background-color: white;
  width: ${Y_AXIS_CHART_WIDTH_PX}px;
`;

export const CellTypeLabelStyle = styled.div`
  margin: auto;
  font: normal 12px sans-serif;
  white-space: pre;
  margin-left: 24px;
  border: none;
  width: 100%;
  color: ${ECHART_AXIS_LABEL_COLOR_HEX};
  text-align: left;
`;

interface TissueHeaderLabelStyleProps extends CommonThemeProps {
  expanded: boolean;
}

export const TissueHeaderLabelStyle = styled.div`
  margin: auto;
  font: normal 12px sans-serif;
  border: none;
  width: 100%;
  color: ${ECHART_AXIS_LABEL_COLOR_HEX};
  text-align: left;
  ${(props: TissueHeaderLabelStyleProps) => {
    return props.expanded && `font-weight: ${fontWeightSemibold(props)}`;
  }}
`;

export const HiddenCellTypeLabelStyle = styled.div`
  /* Overlay invisible, un-intractable element with full name of cell type for ctrl+f page search */
  position: absolute;
  z-index: 1;
  color: rgba(0, 0, 0, 0);
  pointer-events: none;
  user-select: none;
`;

export const CellTypeLabelTooltipStyle = styled.div`
  margin: 0 -10px;
  text-align: center;
`;

export const CellCountLabelStyle = styled.div`
  height: ${HEAT_MAP_BASE_CELL_PX}px;
  background-color: white;
  font: normal normal ${SELECTED_STYLE.fontSize}px ${SELECTED_STYLE.fontFamily};
  white-space: pre;
  border: none;
  color: ${gray500};
  text-align: right;
  padding-top: 3px;
`;
export const FlexRowJustified = styled.div`
  display: flex;
  flex-direction: row;
  justify-content: space-between;
  padding-left: 5px;
  width: 100%;
`;

export const FlexRow = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
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

export const TissueLabel = styled.div`
  white-space: nowrap;
`;
