import styled from "@emotion/styled";
import { fontHeaderL, fontHeaderXl } from "@czi-sds/components";
import { SKINNY_MODE_BREAKPOINT_WIDTH } from "../CellGuideCard/constants";

export const Wrapper = styled.div`
  display: flex;
  flex-direction: column;
  margin-top: 222px;
  width: 400px;
  margin-left: auto;
  margin-right: auto;

  // Mobile landing page styling
  @media (max-width: ${SKINNY_MODE_BREAKPOINT_WIDTH}px) {
    max-width: 100vw;
    padding-left: 24px;
    padding-right: 24px;
    margin-top: 16px;
  }
`;

export const StyledHeader = styled.div`
  ${fontHeaderXl}
  margin: 0 0 16px 8px;

  @media (max-width: ${SKINNY_MODE_BREAKPOINT_WIDTH}px) {
    ${fontHeaderL}
  }
`;
