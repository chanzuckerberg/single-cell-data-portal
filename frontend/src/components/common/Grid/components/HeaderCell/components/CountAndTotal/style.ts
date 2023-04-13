import styled from "@emotion/styled";
import { CommonThemeProps, fontBodyXs, getColors, Tag } from "czifui";

const gray100 = (props: CommonThemeProps) => getColors(props)?.gray[100];
const gray600 = (props: CommonThemeProps) => getColors(props)?.gray[600];

export const CountAndTotal = styled(Tag)`
  &.MuiChip-root {
    background-color: ${gray100};
    margin: 0;

    .MuiChip-label {
      ${fontBodyXs};
      color: ${gray600};
      font-weight: 500;
      letter-spacing: -0.003em;
    }
  }
`;
