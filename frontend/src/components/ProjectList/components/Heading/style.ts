import { fontStyle as fontStyleBase } from "src/components/common/theme";
import { LAYOUT } from "src/components/ProjectList/common/style";
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
`;

export const Info = styled.div`
  ${fontStyle}
  width: ${LAYOUT.INFO}%;
`;

export const QuestionMarkWrapper = styled.span`
  position: absolute;
`;

export const QuestionMark = styled.img`
  position: relative;
  cursor: pointer;
  width: 15px;
  height: 15px;
  top: 3px;
  left: 5px;
`;
