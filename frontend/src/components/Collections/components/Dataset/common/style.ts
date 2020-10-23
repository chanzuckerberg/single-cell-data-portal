import { LAYOUT } from "src/components/Collections/common/layout";
import { BLUE, fontStyle } from "src/components/common/theme";
import styled, { css } from "styled-components";

export const columnStyle = css`
  align-items: center;
  display: flex;
  flex: ${LAYOUT.SMALL_COLUMN};
  justify-content: center;
`;

export const SmallColumn = styled.div`
  ${columnStyle}
  ${fontStyle}
`;

export const buttonStyle = css`
  ${fontStyle}
  text-decoration: none;
  color: ${BLUE};
  background-color: white;
  border: 1px solid ${BLUE};
  border-radius: 3px;
  appearance: none;
  padding: 2px 20px;
  text-align: center;
  box-sizing: border-box;

  &:hover {
    background-color: ${BLUE};
    color: white;
    text-decoration: none;
    cursor: pointer;
  }
`;
