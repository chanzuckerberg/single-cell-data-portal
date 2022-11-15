import { AnchorButton, Button, Classes } from "@blueprintjs/core";
import styled from "@emotion/styled";
import { BLUE, PRIMARY_BLUE, PT_TEXT_COLOR } from "src/components/common/theme";
import { css } from "@emotion/react";

/* Basic button styles. */
const buttonStyle = css`
  &.${Classes.BUTTON} {
    border-radius: 4px;
    box-shadow: none !important; /* required; overrides specificity of bp4 button box shadow rule with important style declaration. */
    font-weight: 500;
    height: 32px; /* overrides BP button height specification defined in global.scss */
    letter-spacing: -0.006em;
    line-height: 20px;

    :focus {
      outline: none;
    }
  }
`;

/* Basic primary button styles. */
const primaryButtonStyle = css`
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

/* Basic styled button. */
const StyledButton = styled(Button)`
  ${buttonStyle}
`;

/* Basic styled anchor button. */
const StyledAnchorButton = styled(AnchorButton)`
  ${buttonStyle}
`;

/* Basic default button with no visual intent color applied to element. */
export const StyledDefaultButton = styled(StyledButton)`
  &.${Classes.BUTTON} {
    background-image: none;
    color: ${PT_TEXT_COLOR};
  }
`;

/* Basic styled primary button. */
export const StyledPrimaryButton = styled(StyledButton)`
  ${primaryButtonStyle}
`;

/* Basic styled outlined button. */
export const StyledOutlineButton = styled(StyledButton)`
  &.${Classes.BUTTON}.${Classes.OUTLINED} {
    /* Primary outlined button. */

    &.${Classes.INTENT_PRIMARY} {
      border-color: ${PRIMARY_BLUE};
      color: ${PRIMARY_BLUE};

      &:hover {
        border-color: ${BLUE.B};
        color: ${BLUE.B};
      }

      &.${Classes.DISABLED} {
        /* maintains BP disabled border color specification; BLUE.B with 20% opacity */
        border-color: rgba(0, 95, 198, 0.2);
        /* maintains BP disabled color specification; BLUE.B with 50% opacity */
        color: rgba(0, 95, 198, 0.5);
      }
    }
  }
`;

/* Basic styled primary minimal button (to be used in conjunction with minimal prop). */
export const StyledPrimaryMinimalButton = styled(StyledButton)`
  &.${Classes.BUTTON}.${Classes.MINIMAL}.${Classes.INTENT_PRIMARY} {
    border-radius: 0; /* overrides basic styled button border radius specification */
    color: ${PRIMARY_BLUE};
    height: unset; /* overrides bp button height specification */
    min-height: unset; /* overrides bp button min height specification */
    padding: 0; /* overrides basic styled button padding specification */

    &:hover {
      background: none;
    }
  }
`;

/* Basic styled primary anchor button. */
export const StyledPrimaryAnchorButton = styled(StyledAnchorButton)`
  ${primaryButtonStyle}
`;
