import { fontStyle } from "src/components/common/theme";
import { LAYOUT } from "src/components/ProjectList/common/layout";
import styled from "styled-components";

export const Wrapper = styled.div`
  width: ${LAYOUT.VIEW}%;
  margin-right: ${LAYOUT.VIEW_MARGIN}%;
  display: flex;
  justify-content: center;
  align-items: center;
`;

export const SelectWrapper = styled.div`
  display: flex;
  position: relative;
`;

export const StyledAnchor = styled.a`
  ${fontStyle}
  text-decoration: none;
  color: #0076dc;
  background-color: white;
  border: 1px solid #0076dc;
  box-shadow: inset 0px -1px 1px rgba(16, 22, 26, 0.2);
  border-radius: 3px;
  appearance: none;
  padding: 0 20px;
  width: 100px;
  height: 30px;
  text-align: center;
`;
