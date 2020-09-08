import { fontStyle } from "src/components/common/theme";
import styled from "styled-components";
import { layout } from "../common/layout";

export const Wrapper = styled.div`
  width: 100%;
  height: 56px;
  padding-left: 15px;
  display: flex;
  align-items: center;
  border-top: 1px solid #dbdcdd;
`;

export const MainWrapper = styled.div`
  ${fontStyle}
  ${layout}
  padding: 0 15px;
  margin: 0 auto;
`;
