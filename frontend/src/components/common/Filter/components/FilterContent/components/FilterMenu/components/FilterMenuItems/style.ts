import styled from "@emotion/styled";
import { CommonThemeProps, getFontWeights } from "@czi-sds/components";

interface Props extends CommonThemeProps {
  selected: boolean;
}

const semibold = (props: CommonThemeProps) => getFontWeights(props)?.semibold; // font-weight 600.

export const ListItem = styled("li")<Props>`
  .MuiListItemText-root {
    span:first-of-type {
      font-weight: ${({ selected }) => (selected ? semibold : undefined)};
    }
  }
`;
