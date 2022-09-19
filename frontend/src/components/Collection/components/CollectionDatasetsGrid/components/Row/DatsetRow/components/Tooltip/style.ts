import styled from "@emotion/styled";
import { fontStyle } from "src/components/common/theme";

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
