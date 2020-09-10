import styled from "styled-components";
import { layout } from "../common/layout";

export const Wrapper = styled.div`
  display: flex;
  position: sticky;
  box-shadow: 0px 2px 4px rgba(0, 0, 0, 0.1);
  top: 0;
  background-color: #f5f5f5;
  height: 60px;
  width: 100%;
  align-items: center;
  /* As required in mock */
  margin-bottom: 32px;
  z-index: 1;
`;

export const MainWrapper = styled.div`
  ${layout}
  padding-left: 22px;
  margin: 0 auto;
`;
