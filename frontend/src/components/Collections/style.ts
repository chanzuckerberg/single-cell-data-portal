import styled from "styled-components";
import { PT_GRID_SIZE_PX } from "../common/theme";

export const TitleWrapper = styled.div`
  display: flex;
  flex-direction: row;
  justify-content: space-between;
  align-items: flex-end;
`;

export const TitleAndDescription = styled.div`
  display: flex;
  flex-direction: column;
  & > h1 {
    margin-bottom: ${PT_GRID_SIZE_PX}px;
  }
  & > p {
    margin-bottom: 0;
  }
`;
