import { css } from "@emotion/react";
import styled from "@emotion/styled";
import { LAYOUT } from "src/components/Collections/common/layout";
import { fontStyle, OLD_BLUE } from "src/components/common/theme";

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
  color: ${OLD_BLUE};
  background-color: white;
  border: 1px solid ${OLD_BLUE};
  border-radius: 3px;
  appearance: none;
  padding: 2px 20px;
  text-align: center;
  box-sizing: border-box;

  &:hover {
    background-color: ${OLD_BLUE};
    color: white;
    text-decoration: none;
    cursor: pointer;
  }
`;

export const disabledButtonStyle = css`
  color: rgb(0, 118, 220, 0.5);
  border: 1px solid rgb(0, 118, 220, 0.2);

  &:hover {
    cursor: not-allowed;
    background-color: white;
  }
`;
