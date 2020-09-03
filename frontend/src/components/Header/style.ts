import { LAYOUT } from "src/components/common/layout";
import styled from "styled-components";

export const Wrapper = styled.div`
  display: flex;
  position: sticky;
  top: 0;
  background-color: #f5f5f5;
  height: 70px;
  width: 100%;
  align-items: center;
  /* As required in mock */
  margin-bottom: 42px;
`;

export const MainWrapper = styled.div`
  max-width: 1400px;
  min-width: 800px;
  width: 100%;
  padding-left: 22px;
  margin: 0 auto;
`;

