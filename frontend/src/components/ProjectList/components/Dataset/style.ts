import { fontStyle } from "src/components/common/theme";
import { LAYOUT } from "src/components/ProjectList/common/layout";
import styled from "styled-components";

export const Wrapper = styled.div`
  :nth-child(even) {
    background: rgba(191, 204, 214, 0.2);
  }
  :last-of-type {
    border-bottom: none;
  }
  display: flex;
  min-height: 40px;
  border-bottom: 1px solid #e5e5e5;
`;

export const Name = styled.div`
  ${fontStyle}
  width: ${LAYOUT.NAME}%;
  margin-right: ${LAYOUT.NAME_MARGIN}%;
  display: flex;
  align-items: center;
  padding-left: 10px;
`;
