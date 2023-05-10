import styled from "@emotion/styled";
import { fontBodyS, getColors } from "czifui";

export const WmgLink = styled.a`
  ${fontBodyS}
  font-weight: 500;
  ${(props) => {
    const colors = getColors(props);
    return `color: ${colors?.primary[400]}`;
  }}
`;
