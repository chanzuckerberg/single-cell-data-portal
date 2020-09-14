import { BLUE, fontStyle } from "src/components/common/theme";
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
  color: ${BLUE};
  background-color: white;
  border: 1px solid ${BLUE};
  border-radius: 3px;
  appearance: none;
  padding: 2px 20px;
  text-align: center;
  box-sizing: border-box;

  &:hover {
    background-color: ${BLUE};
    color: white;
  }
`;
