import { fontStyle } from "src/components/common/theme";
import { LAYOUT } from "src/components/ProjectList/common/layout";
import styled from "styled-components";

export const Wrapper = styled.div`
  width: ${LAYOUT.INFO}%;
  align-items: center;
  line-height: 1.2;
  padding: 12px 0;
  display: flex;
`;

export const StyledAnchor = styled.a`
  ${fontStyle}
  text-decoration: none;
  color: #318ad8;
  display: contents;

  &:hover {
    text-decoration: underline;
  }
`;
