import { Alert, Classes } from "@blueprintjs/core";
import { DARK_GRAY } from "src/components/common/theme";
import styled from "styled-components";

export const StyledAlert = styled(Alert)`
  .${Classes.BUTTON} {
    :not(.${Classes.INTENT_DANGER}) {
      background: none;
      box-shadow: none !important;
      color: ${DARK_GRAY.A};

      &:hover {
        background: rgba(167, 182, 194, 0.3);
      }
    }
  }
`;
