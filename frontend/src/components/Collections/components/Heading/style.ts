import { fontStyle as fontStyleBase } from "src/components/common/theme";
import styled, { css } from "styled-components";
import { LAYOUT } from "../../common/layout";

const fontStyle = css`
  ${fontStyleBase}
  line-height: 21px;
  font-weight: 600;
`;

export const Wrapper = styled.div`
  display: flex;
  padding-bottom: 9px;
  border-bottom: 1px solid #e5e5e5;
`;

export const Name = styled.div`
  ${fontStyle}
  flex: ${LAYOUT.NAME};
  /* (thuang): width is ignored when flex is used */
  min-width: 620px;
  max-width: 620px;
`;

export const QuestionMark = styled.img`
  position: relative;
  cursor: pointer;
  width: 15px;
  height: 15px;
  top: 3px;
  left: 5px;
`;
