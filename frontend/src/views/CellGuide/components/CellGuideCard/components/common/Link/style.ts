import styled from "@emotion/styled";
import { fontBodyS } from "@czi-sds/components";
import { primary400, primary500 } from "src/common/theme";

export const StyledLink = styled.a`
  ${fontBodyS}

  color: ${primary400};

  &:hover {
    text-decoration: none;
    color: ${primary500};
  }
`;
