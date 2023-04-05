import styled from "@emotion/styled";
import { CommonThemeProps, fontBodyXs, getColors, Tag } from "czifui";

const gray100 = (props: CommonThemeProps) => getColors(props)?.gray[100];

export const CountAndTotal = styled(Tag)`
  &.MuiChip-root {
    background-color: ${gray100};
    cursor: inherit;
    margin: 0;

    .MuiChip-label {
      ${fontBodyXs};
      color: #000000;
      font-weight: 500;
      letter-spacing: -0.003em;
    }
  }
`;
