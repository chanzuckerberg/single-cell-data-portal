import { fontBodyS } from "@czi-sds/components";
import styled from "@emotion/styled";
import { primary400 } from "src/common/theme";

export const MarkerGeneTableWrapper = styled.section`
  margin-top: -40px;
`;

export const StyledLink = styled.a`
  ${fontBodyS}
  color: ${primary400};
  &:hover {
    text-decoration: none !important;
  }
`;
