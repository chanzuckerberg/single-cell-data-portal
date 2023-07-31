import styled from "@emotion/styled";
import { CommonThemeProps } from "@czi-sds/components";
import { fontWeightSemibold } from "src/common/theme";

interface Props extends CommonThemeProps {
  selected: boolean;
}

export const ListItem = styled("li")<Props>`
  .MuiListItemText-root {
    span:first-of-type {
      font-weight: ${(props) =>
        props.selected ? fontWeightSemibold(props) : undefined};
    }
  }
`;
