import styled from "@emotion/styled";
import { fontBodyXxxs, getColors } from "@czi-sds/components";

export const StyledLegendText = styled.text`
  ${fontBodyXxxs}

  alignment-baseline: middle;

  ${(props) => {
    const colors = getColors(props);
    return `color: ${colors?.gray[500]}`;
  }}
`;
