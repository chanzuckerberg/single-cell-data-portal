import { LAYOUT } from "src/components/ProjectList/common/style";
import styled from "styled-components";

export const Wrapper = styled.div`
  width: ${LAYOUT.VIEW}%;
  margin-right: ${LAYOUT.VIEW_MARGIN}%;
  display: flex;
  justify-content: center;
  align-items: center;
`;

export const StyledAnchor = styled.a`
  text-decoration: none;
  font-size: 12px;
  line-height: 18px;
  color: #318ad8;
`;
