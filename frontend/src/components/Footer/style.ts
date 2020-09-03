import { fontStyle } from "src/components/common/theme";
import { LAYOUT } from "src/components/common/layout";
import styled from "styled-components";

export const Wrapper = styled.div`
  width: 100%;
  margin-right: ${LAYOUT.NAME_MARGIN}%;
  height: 56px;
  padding-left: 15px;
  display: flex;
  align-items: center;
  border-top: 1px solid #DBDCDD;
`;

export const MainWrapper = styled.div`
  ${fontStyle}
  max-width: 1400px;
  min-width: 800px;
  width: 100%;
  padding: 0 15px;
  margin: 0 auto;
`;
