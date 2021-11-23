import styled from "styled-components";
import { layout } from "../common/layout";

export const Wrapper = styled.div`
  display: flex;
  flex-direction: column;
  min-height: 100vh;
`;

export const MainWrapper = styled.main`
  ${layout}
  min-height: 83vh;
  padding: 15px 25px;
  margin: 0 auto 32px auto;
  display: flex;
  flex: 1;
`;
