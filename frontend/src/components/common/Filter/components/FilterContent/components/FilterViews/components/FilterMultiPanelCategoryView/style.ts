import { Classes } from "@blueprintjs/core";
import styled from "@emotion/styled";
import FilterSearch from "src/components/common/Filter/components/FilterSearch";
import { spacesXs } from "src/common/theme";

export const Views = styled("div")`
  display: flex;
  flex-direction: column;
  padding: ${spacesXs}px;
`;

export const Search = styled(FilterSearch)`
  padding: ${spacesXs}px;

  &.${Classes.INPUT_GROUP} {
    .${Classes.ICON} {
      left: ${spacesXs}px;
    }

    .${Classes.INPUT_ACTION} {
      right: 14px;
    }
  }
`;
