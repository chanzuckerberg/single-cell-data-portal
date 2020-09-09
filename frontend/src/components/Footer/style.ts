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
  padding-left: 15px;
  margin: 0 auto;
`;

export const StyledAnchor = styled.a`
  ${fontStyle}
  text-decoration: none;
  color: #0076dc;
  appearance: none;
  padding-left: 15px;
  width: 100px;
  height: 30px;
  text-align: center;
`;
