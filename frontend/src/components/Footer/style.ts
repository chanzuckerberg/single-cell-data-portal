import styled from "@emotion/styled";
import { fontStyle, OLD_BLUE } from "src/components/common/theme";
import { layout } from "../common/layout";

export const FOOTER_HEIGHT_PX = 56;

export const Wrapper = styled.div`
  width: 100%;
  height: ${FOOTER_HEIGHT_PX}px;
  padding-left: 15px;
  display: flex;
  flex: none;
  align-items: center;
  border-top: 1px solid #dbdcdd;
  z-index: 1;
  background: white;
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
  color: ${OLD_BLUE};
  appearance: none;
  padding-left: 15px;
  width: 100px;
  height: 30px;
  text-align: center;
`;
