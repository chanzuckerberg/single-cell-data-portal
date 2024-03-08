import styled from "@emotion/styled";
import Link from "next/link";
import { primary400, primary500 } from "src/common/theme";

export const StyledAnchor = styled(Link)`
  color: ${primary400};
  display: block;

  &:focus {
    outline: none;
  }

  &:hover {
    background: transparent;
    color: ${primary500};
    text-decoration: none;
  }
`;
