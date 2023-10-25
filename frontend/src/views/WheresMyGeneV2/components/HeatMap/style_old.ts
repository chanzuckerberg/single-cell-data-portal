import styled from "@emotion/styled";
import { Y_AXIS_CHART_WIDTH_PX } from "./utils";
import { LIGHT_GRAY } from "src/components/common/theme";
import { LEGEND_HEIGHT_PX } from "../../../WheresMyGene/components/InfoPanel/components/Legend/style";
import { HEADER_HEIGHT_PX } from "src/components/Header/style";
import {
  CONTENT_WRAPPER_LEFT_RIGHT_PADDING_PX,
  CONTENT_WRAPPER_TOP_BOTTOM_PADDING_PX,
} from "src/components/Layout/style";
import { LEGEND_MARGIN_BOTTOM_PX } from "../../../WheresMyGene/style";
import { X_AXIS_CHART_HEIGHT_PX } from "../../common/constants";
import { Autocomplete } from "@mui/material";
import { TagFilter } from "@czi-sds/components";
import { spacesS } from "src/common/theme";

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

export const StyledAutocomplete = styled(Autocomplete)`
  width: 258px;

  .MuiInputBase-root {
    padding-right: 8px !important;
  }
  & .MuiAutocomplete-inputRoot[class*="MuiInput-root"] {
    padding: 0px;
  }
  &
    .MuiAutocomplete-inputRoot[class*="MuiOutlinedInput-root"]
    .MuiAutocomplete-input {
    padding: 0px;
  }
  & .MuiInputLabel-root {
    margin-top: -8px;
  }
  & .MuiInputLabel-root.MuiInputLabel-shrink {
    margin-top: 0px;
  }
`;

export const StyledTag = styled(TagFilter)`
  max-width: 258px;
`;
// Copied to frontend/src/views/WheresMyGeneV2/components/HeatMap/style.ts
interface TopLeftCornerMaskProps {
  height: number;
}

// Copied to frontend/src/views/WheresMyGeneV2/components/HeatMap/style.ts
export const TopLeftCornerMask = styled.div<TopLeftCornerMaskProps>`
  position: absolute;
  background-color: white;
  z-index: 3;
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

/**
 * (thuang): Instead of using the full width, we only want enough space for the
 * filter box + the widest scrollbar across different browsers.
 */
/**
 * Copied to frontend/src/views/WheresMyGene/components/HeatMap/utils.ts for
 * frontend/src/views/WheresMyGeneV2/components/HeatMap/style.ts to import
 */
const CELL_TYPE_FILTER_WIDTH_PX = 300;

// Copied to frontend/src/views/WheresMyGeneV2/components/HeatMap/style.ts
export const CellTypeFilterContainer = styled.div`
  height: 100%;
  width: ${CELL_TYPE_FILTER_WIDTH_PX}px;
`;

const CELL_TYPE_SEARCH_BOX_HEIGHT_PX = 37;

export const CellTypeTagContainer = styled.div`
  overflow-y: auto;
  height: calc(100% - ${CELL_TYPE_SEARCH_BOX_HEIGHT_PX}px);
  padding: ${spacesS}px;
`;

interface YAxisWrapperProps {
  top: number;
}

export const YAxisWrapper = styled.div<YAxisWrapperProps>`
  width: ${Y_AXIS_CHART_WIDTH_PX}px;
  position: sticky;
  top: ${(props) => props.top ?? X_AXIS_CHART_HEIGHT_PX}px;
  left: 0;
  z-index: 1;
  padding-top: 5px;
  /* Somehow Firefox requires this to scroll */
  overflow: hidden;
`;

interface XAxisMaskProps {
  height: number;
}

export const XAxisMask = styled.div<XAxisMaskProps>`
  width: ${Y_AXIS_CHART_WIDTH_PX + CHART_PADDING_PX}px;
  height: ${(props) => props.height}px;
`;

// Copied to frontend/src/views/WheresMyGeneV2/components/HeatMap/style.ts
export const XAxisWrapper = styled.div`
  display: flex;
  background-color: white;
  flex-direction: row;
  position: sticky;
  top: 0;
  width: 100%;
  z-index: 2;
`;

interface ChartWrapperProps {
  top: number;
}

export const ChartWrapper = styled.div<ChartWrapperProps>`
  position: absolute;
  padding-left: ${CHART_PADDING_PX}px;
  padding-right: ${CHART_PADDING_PX}px;
  padding-top: 5px;
  left: ${Y_AXIS_CHART_WIDTH_PX}px;
  top: ${(props) => props.top}px;
`;
