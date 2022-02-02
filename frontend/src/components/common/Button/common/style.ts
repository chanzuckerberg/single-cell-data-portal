import { Button, Classes } from "@blueprintjs/core";
import { BLUE, PRIMARY_BLUE } from "src/components/common/theme";
import styled from "styled-components";

/* Basic styled button. */
const StyledButton = styled(Button)`
  &.${Classes.BUTTON} {
    border-radius: 3px;
    box-shadow: none !important; /* required; overrides specificity of bp3 button box shadow rule with important style declaration. */
    font-weight: 500;
    letter-spacing: -0.1px;
    line-height: 18px;
  }
`;

/* Basic styled primary button. */
export const StyledPrimaryButton = styled(StyledButton)`
  &.${Classes.BUTTON}.${Classes.INTENT_PRIMARY} {
    background-color: ${PRIMARY_BLUE};
    &:hover {
      /* maintains BP hover background color specification */
      background-color: ${BLUE.B};
    }
    &.${Classes.DISABLED} {
      /* maintains BP disabled background color specification; BLUE.C with 50% opacity */
      background-color: rgba(14, 125, 236, 0.5);
    }
  }
`;

/* Basic styled primary outlined button. */
export const StyledOutlineButton = styled(StyledButton)`
  &.${Classes.BUTTON}.${Classes.OUTLINED}.${Classes.INTENT_PRIMARY} {
    border-color: ${PRIMARY_BLUE};
    color: ${PRIMARY_BLUE};
    &:hover {
      /* maintains BP hover color specification */
      color: ${BLUE.B};
    }
    &.${Classes.DISABLED} {
      /* maintains BP disabled border color specification; BLUE.B with 20% opacity */
      border-color: rgba(0, 95, 198, 0.2);
      /* maintains BP disabled color specification; BLUE.B with 50% opacity */
      color: rgba(0, 95, 198, 0.5);
    }
  }
`;
