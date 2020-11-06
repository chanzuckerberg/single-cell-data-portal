import styled from "styled-components";
import { layout } from "../common/layout";
import { PT_GRID_SIZE_PX } from "../common/theme";

export const Wrapper = styled.div`
  display: flex;
  position: sticky;
  box-shadow: 0px 0px 0px rgba(16, 22, 26, 0.1),
    0px 1px 1px rgba(16, 22, 26, 0.2);
  top: 0;
  background-color: white;
  height: 60px;
  width: 100%;
  align-items: center;
  /* As required in mock */
  margin-bottom: 32px;
  z-index: 1;
`;

export const MainWrapper = styled.div`
  ${layout}
  display: flex;
  justify-items: space-between;
  padding-left: 22px;
  margin: 0 auto;
  justify-content: space-between;
  align-items: center;
`;

export const ButtonWrapper = styled.div`
  a:not(:last-child) {
    margin-right: ${2 * PT_GRID_SIZE_PX}px;
  }

  margin-right: ${3 * PT_GRID_SIZE_PX}px;
`;
