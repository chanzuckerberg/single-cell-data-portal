import styled from "@emotion/styled";
import {
  CommonThemeProps,
  fontBodyXs,
  getColors,
  getFontWeights,
} from "czifui";

const fontWeight600 = (props: CommonThemeProps) =>
  getFontWeights(props)?.semibold;
const primary400 = (props: CommonThemeProps) => getColors(props)?.primary[400];

export const SelectedTags = styled.span`
  display: flex;
  flex-wrap: wrap;
  gap: 8px; /* TODO(cc) confirm */
  min-width: 0; /* facilitates ellipsis on tag should it be required; flex default for min width is "auto" */

  .MuiChip-root {
    margin: 0;

    &:hover {
      background-color: ${primary400};
    }

    .MuiChip-label {
      ${fontBodyXs};
      font-weight: ${fontWeight600};
      letter-spacing: -0.003em;
      white-space: normal;
    }

    .MuiSvgIcon-root {
      height: 16px;
      width: 16px;
    }
`;
