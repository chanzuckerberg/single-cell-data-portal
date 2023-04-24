import styled from "@emotion/styled";
import { CommonThemeProps, getColors } from "czifui";

const primary400 = (props: CommonThemeProps) => getColors(props)?.primary[400];
const primary500 = (props: CommonThemeProps) => getColors(props)?.primary[500];

export const StyledAnchor = styled.a`
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
