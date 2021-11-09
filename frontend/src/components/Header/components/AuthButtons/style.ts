import { BLUE, PT_GRID_SIZE_PX } from "src/components/common/theme";
import styled from "styled-components";

export const Initial = styled.div`
  align-items: center;
  background-color: ${BLUE.F};
  border-radius: ${3 * PT_GRID_SIZE_PX}px;
  display: flex;
  font-weight: 600;
  height: ${3 * PT_GRID_SIZE_PX}px;
  justify-content: center;
  text-transform: capitalize;
  width: ${3 * PT_GRID_SIZE_PX}px;
`;
