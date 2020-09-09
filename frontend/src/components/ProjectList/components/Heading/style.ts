import { fontStyle as fontStyleBase } from "src/components/common/theme";
import { LAYOUT } from "src/components/ProjectList/common/layout";
import styled, { css } from "styled-components";

const fontStyle = css`
  ${fontStyleBase}
  line-height: 21px;
  font-weight: 600;
`;

export const Wrapper = styled.div`
  display: flex;
  padding-bottom: 9px;
  border-bottom: 1px solid black;
  box-shadow: 0px 2px 4px rgba(0, 0, 0, 0.1)
`;

export const Name = styled.div`
  ${fontStyle}
  margin-right: ${LAYOUT.NAME_MARGIN}%;
  width: ${LAYOUT.NAME}%;
`;

export const View = styled.div`
  ${fontStyle}
  margin-right: ${LAYOUT.VIEW_MARGIN}%;
  width: ${LAYOUT.VIEW}%;
  display: flex;
  align-items: center;
  justify-content: center;
`;

export const Info = styled.div`
  ${fontStyle}
  width: ${LAYOUT.INFO}%;
`;

export const QuestionMark = styled.img`
  position: relative;
  cursor: pointer;
  width: 15px;
  height: 15px;
  top: 3px;
  left: 5px;
`;
