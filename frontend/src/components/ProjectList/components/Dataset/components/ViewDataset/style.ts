import { BLUE, fontStyle } from "src/components/common/theme";
import styled from "styled-components";

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
