import { PRIMARY_BLUE } from "src/components/common/theme";
import styled from "styled-components";

export const StyledAnchor = styled.a`
  color: inherit;

  &:focus {
    outline: none;
  }

  &:hover {
    background: transparent;
    color: ${PRIMARY_BLUE};
    text-decoration: none;
  }
`;
