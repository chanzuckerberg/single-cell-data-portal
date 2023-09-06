import styled from "@emotion/styled";
import {
  CommonThemeProps,
  fontHeaderL,
  fontHeaderXl,
} from "@czi-sds/components";

interface WrapperProps extends CommonThemeProps {
  searchBarOpen: boolean;
}

export const Wrapper = styled.div<WrapperProps>`
  display: flex;
  flex-direction: column;
  margin-top: 222px;
  width: 400px;
  margin-left: auto;
  margin-right: auto;

  // Mobile landing page styling
  @media (max-width: 768px) {
    max-width: 100vw;
    padding-left: 24px;
    padding-right: 24px;
    justify-content: ${({ searchBarOpen }) =>
      searchBarOpen ? "flex-start" : "center"};
    height: 100vh;
    margin-top: ${({ searchBarOpen }) => (searchBarOpen ? "0" : "-10vh")};
  }
`;

export const StyledHeader = styled.div`
  ${fontHeaderXl}
  margin: 0 0 16px 8px;

  @media (max-width: 768px) {
    ${fontHeaderL}
  }
`;
