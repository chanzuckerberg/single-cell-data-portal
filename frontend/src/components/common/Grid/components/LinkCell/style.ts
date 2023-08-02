import styled from "@emotion/styled";
import { primary400, primary500 } from "src/common/theme";

export const StyledAnchor = styled.a`
  color: ${primary400};
  display: block;

  &:focus {
    outline: none;
  }

  &:hover {
    background: transparent;
    color: ${primary500};
    text-decoration: none;
  }
`;
