import { fontStyle } from "src/components/common/theme";
import styled from "styled-components";

export const StyledButton = styled.button`
  ${fontStyle}
  text-decoration: none;
  color: #0076dc;
  background-color: white;
  border: 1px solid #0076dc;
  border-radius: 3px;
  appearance: none;
  padding: 2px 20px;
  text-align: center;
  box-sizing: border-box;

  &:hover {
    background-color: #0076dc;
    color: white;
  }
`;
