import { Alert, Classes } from "@blueprintjs/core";
import { BLUE, RED } from "src/components/common/theme";
import styled from "styled-components";

export const StyledAlert = styled(Alert)`
  .${Classes.BUTTON} {
    background: ${RED.C};
    box-shadow: none !important;

    :not(.${Classes.INTENT_DANGER}) {
      color: ${BLUE.C};
      background: none;
    }
  }
`;
