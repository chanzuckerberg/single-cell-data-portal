import styled from "@emotion/styled";
import { fontBodyXxs } from "@czi-sds/components";
import { gray500 } from "src/common/theme";
import { Button } from "@czi-sds/components";
import {
  CONTENT_WRAPPER_LEFT_RIGHT_PADDING_PX,
  CONTENT_WRAPPER_TOP_BOTTOM_PADDING_PX,
} from "src/components/Layout/style";
import { LEGEND_HEIGHT_PX } from "src/views/WheresMyGeneV2/components/InfoPanel/components/Legend/style";
import { Y_AXIS_CHART_WIDTH_PX } from "src/views/WheresMyGeneV2/components/HeatMap/utils";
import {
  GENE_SEARCH_BAR_HEIGHT_PX,
  GENE_SEARCH_LEFT_OFFSET_PX,
} from "src/views/WheresMyGeneV2/common/constants";
import { HEADER_HEIGHT_PX } from "src/components/Header/style";

interface ContainerProps {
  sidebarWidth: number;
}

export const Container = styled.div`
  height: ${GENE_SEARCH_BAR_HEIGHT_PX}px;
  width: fit-content;
  position: fixed;
  top: ${HEADER_HEIGHT_PX +
  CONTENT_WRAPPER_TOP_BOTTOM_PADDING_PX +
  LEGEND_HEIGHT_PX}px;
  /**
   * (thuang): Dynamically calculate the left offset based on the sidebar
   * expand/collapse width.
   */
  left: ${({ sidebarWidth }: ContainerProps) =>
    sidebarWidth +
    CONTENT_WRAPPER_LEFT_RIGHT_PADDING_PX +
    Y_AXIS_CHART_WIDTH_PX +
    GENE_SEARCH_LEFT_OFFSET_PX}px;
`;

export const AutocompleteWrapper = styled.div`
  width: 240px;
`;

export const ActionWrapper = styled.div`
  display: flex;
  gap: 16px;
`;

export const Label = styled.label`
  ${fontBodyXxs}

  color: ${gray500};
`;

export const LoadingIndicatorWrapper = styled.div`
  display: flex;
  align-items: center;
`;

export const StyledButtonWrapper = styled.div`
  align-self: center;
`;

export const StyledClearButton = styled(Button)`
  white-space: nowrap;
  font-weight: 500;
`;
