import { fontStyle } from "src/components/common/theme";
import { LAYOUT } from "src/components/ProjectList/common/style";
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
  color: #318ad8;
  background-color: white;
  border: 1px solid #318ad8;
  border-radius: 5px;
  appearance: none;
  padding: 0 20px;
  width: 180px;
  height: 30px;
  text-align: center;
`;
