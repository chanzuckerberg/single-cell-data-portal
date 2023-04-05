import styled from "@emotion/styled";
import { CommonThemeProps, getColors } from "czifui";

const primary400 = (props: CommonThemeProps) => getColors(props)?.primary[400];

export const StyledAnchor = styled.a`
  color: ${primary400};
  font-weight: 500;

  &:focus {
    outline: none;
  }

  &:hover {
    background: transparent;
    color: ${primary400};
    text-decoration: none;
  }
`;
