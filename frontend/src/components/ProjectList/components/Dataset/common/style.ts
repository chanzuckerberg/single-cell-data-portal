import { fontStyle } from "src/components/common/theme";
import { LAYOUT } from "src/components/ProjectList/common/layout";
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
