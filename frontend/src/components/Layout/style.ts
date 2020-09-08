import styled from "styled-components";
import { layout } from "../common/layout";

export const Wrapper = styled.div`
  display: flex;
  flex-direction: column;
`;

export const MainWrapper = styled.div`
  ${layout}
  min-height: 83vh;
  padding: 15px 25px;
  margin: 0 auto;
`;
