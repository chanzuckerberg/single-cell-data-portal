import { PRIMARY_BLUE, PRIMARY_BLUE_500 } from "src/components/common/theme";
import styled from "styled-components";

export const StyledAnchor = styled.a`
  color: ${PRIMARY_BLUE};

  &:focus {
    outline: none;
  }

  &:hover {
    background: transparent;
    color: ${PRIMARY_BLUE_500};
    text-decoration: none;
  }
`;
