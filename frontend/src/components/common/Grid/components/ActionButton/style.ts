import { Button, Classes } from "@blueprintjs/core";
import { PRIMARY_BLUE } from "src/components/common/theme";
import styled from "styled-components";

export const ActionButton = styled(Button)`
  &.${Classes.BUTTON}.${Classes.MINIMAL} {
    height: 24px;
    padding: 4px;
    width: 24px;

    svg {
      fill: ${PRIMARY_BLUE};
    }

    &.${Classes.DISABLED}, &:disabled {
      color: rgba(92, 112, 128, 0.6);

      svg {
        fill: rgba(92, 112, 128, 0.6);
      }
    }

    &:focus {
      outline: none;
    }
  }
`;
