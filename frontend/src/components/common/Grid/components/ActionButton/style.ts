import { Button, Classes } from "@blueprintjs/core";
import styled from "styled-components";

export const ActionButton = styled(Button)`
  &.${Classes.BUTTON}.${Classes.MINIMAL} {
    height: 32px;
    padding: 8px;
    width: 32px;

    &:focus {
      outline: none;
    }
  }
`;
