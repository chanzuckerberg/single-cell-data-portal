import styled from "@emotion/styled";
import { CommonThemeProps, fontBodyXs, getColors, getSpaces } from "czifui";

const primary400 = (props: CommonThemeProps) => getColors(props)?.primary[400];
const spacesS = (props: CommonThemeProps) => getSpaces(props)?.s;

export const SelectedTags = styled.span`
  display: flex;
  flex-wrap: wrap;
  gap: ${spacesS}px;
  min-width: 0; /* facilitates ellipsis on tag should it be required; flex default for min width is "auto" */

  .MuiChip-root {
    margin: 0;

    &:hover,
    &:active {
      background-color: ${primary400};
    }

    &:active {
      box-shadow: none;
    }

    &:focus:hover,
    &:focus:active {
      background-color: ${primary400};
    }

    .MuiChip-label {
      ${fontBodyXs}
      font-weight: 500;
      letter-spacing: -0.003em;
      white-space: normal;
    }

    .MuiSvgIcon-root {
      height: 10px;
      width: 10px;
    }
  }
`;
