import { BLUE, PT_GRID_SIZE_PX } from "src/components/common/theme";
import styled from "styled-components";

export const ButtonWrapper = styled.div`
  a:not(:last-child) {
    margin-right: ${2 * PT_GRID_SIZE_PX}px;
  }

  margin-right: ${3 * PT_GRID_SIZE_PX}px;
`;

export const Initial = styled.div`
  width: ${3 * PT_GRID_SIZE_PX}px;
  height: ${3 * PT_GRID_SIZE_PX}px;
  background-color: ${BLUE.F};
  border-radius: ${3 * PT_GRID_SIZE_PX}px;
  text-transform: capitalize;
  font-weight: 600;
  display: flex;
  align-items: center;
  justify-content: center;
`;
