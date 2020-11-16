import { PT_GRID_SIZE_PX } from "src/components/common/theme";
import styled from "styled-components";

export const ButtonWrapper = styled.div`
  a:not(:last-child) {
    margin-right: ${2 * PT_GRID_SIZE_PX}px;
  }

  margin-right: ${3 * PT_GRID_SIZE_PX}px;
`;

export const Scientist = styled.span`
  font-size: ${2 * PT_GRID_SIZE_PX}px;
`;
