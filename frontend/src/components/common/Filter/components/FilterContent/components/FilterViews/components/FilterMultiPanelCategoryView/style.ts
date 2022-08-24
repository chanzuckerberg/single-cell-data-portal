import { Classes } from "@blueprintjs/core";
import styled from "@emotion/styled";
import FilterSearch from "src/components/common/Filter/components/FilterSearch";

export const Views = styled("div")`
  padding: 6px;
`;

export const Search = styled(FilterSearch)`
  padding: 6px;

  &.${Classes.INPUT_GROUP} {
    .${Classes.ICON} {
      left: 6px;
    }
  }
`;
