import { PT_GRID_SIZE_PX } from "src/components/common/theme";
import styled from "styled-components";

export const Wrapper = styled.div`
  grid-column: 6 / span 3;
  justify-self: end;
  display: flex;
  flex-direction: row;
  & > :not(:last-child) {
    margin-right: ${PT_GRID_SIZE_PX * 2}px;
  }
`;
