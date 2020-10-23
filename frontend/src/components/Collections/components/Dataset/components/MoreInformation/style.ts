import { fontStyle } from "src/components/common/theme";
import styled from "styled-components";

export const StyledAnchor = styled.a`
  ${fontStyle}
  text-decoration: none;
  color: #318ad8;
  display: contents;

  &:hover {
    text-decoration: underline;
  }
`;

export const Wrapper = styled.div`
  display: flex;
  justify-content: center;
  align-items: center;
  flex-direction: column;
`;
