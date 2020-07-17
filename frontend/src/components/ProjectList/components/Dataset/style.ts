import { fontStyle } from "src/components/common/theme";
import { LAYOUT } from "src/components/ProjectList/common/style";
import styled from "styled-components";

export const Wrapper = styled.div`
  display: flex;
  min-height: 60px;
  border-bottom: 1px solid #e5e5e5;
`;

export const Name = styled.div`
  ${fontStyle}
  width: ${LAYOUT.NAME}%;
  margin-right: ${LAYOUT.NAME_MARGIN}%;
  display: flex;
  align-items: center;
`;
