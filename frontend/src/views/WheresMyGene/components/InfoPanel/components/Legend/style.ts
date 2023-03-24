import styled from "@emotion/styled";

export const LEGEND_HEIGHT_PX = 52;
export const LEGEND_MARGIN_BOTTOM_PX = 20;

export const LegendWrapper = styled.div`
  display: flex;
  width: 700px;
  justify-content: flex-end;
  height: ${LEGEND_HEIGHT_PX}px;
  margin-bottom: ${LEGEND_MARGIN_BOTTOM_PX}px;
`;
