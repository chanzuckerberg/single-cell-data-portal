import { PT_GRID_SIZE_PX } from "src/components/common/theme";
import styled from "styled-components";

export const StyledDiv = styled.div`
  display: flex;
  flex-direction: row;
  & > :not(:last-child) {
    margin-right: ${PT_GRID_SIZE_PX}px;
  }
`;
