import styled from "@emotion/styled";
import { fontBodyS, getColors } from "czifui";

export const StyledLink = styled.a`
  ${fontBodyS}

  ${(props) => {
    const colors = getColors(props);
    return `
      color: ${colors?.primary[400]};
    `;
  }}

  &:hover {
    text-decoration: none;
    ${(props) => {
      const colors = getColors(props);
      return `
        color: ${colors?.primary[500]};
      `;
    }}
  }
`;
