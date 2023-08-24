import styled from "@emotion/styled";
import { fontHeaderXxl } from "@czi-sds/components";

export const Wrapper = styled.div`
  display: flex;
  flex-direction: column;
  margin-top: 222px;
  width: 400px;
  margin-left: auto;
  margin-right: auto;

  // Adds padding to cell guide landing page
  @media (max-width: 768px) {
    width: 90vw;
  }
`;

export const StyledHeader = styled.div`
  ${fontHeaderXxl}
  margin-bottom: 24px;
`;
