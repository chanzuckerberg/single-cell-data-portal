import { fontStyle } from "src/components/common/theme";
import { LAYOUT } from "src/components/ProjectList/common/layout";
import styled from "styled-components";

export const Wrapper = styled.div`
  width: ${LAYOUT.INFO}%;
  align-items: center;
`;

export const StyledAnchor = styled.a`
  ${fontStyle}
  text-decoration: none;
  line-height: 18px;
  color: #318ad8;
`;
