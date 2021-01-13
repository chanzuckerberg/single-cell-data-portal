import { fontStyle } from "src/components/common/theme";
import styled from "styled-components";

export const StyledAnchor = styled.a`
  ${fontStyle}
  text-decoration: none;
  color: white;
  display: contents;

  &:hover {
    color: white;
    text-decoration: underline;
  }
`;
