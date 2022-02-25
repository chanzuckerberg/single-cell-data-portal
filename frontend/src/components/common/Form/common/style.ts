import { Classes, Label } from "@blueprintjs/core";
import {
  GRAY,
  LIGHT_GRAY,
  PRIMARY_BLUE,
  PT_TEXT_COLOR,
  RED,
} from "src/components/common/theme";
import styled, { css } from "styled-components";

export const formField = css`
  border-radius: 3px;
  color: ${PT_TEXT_COLOR};
  height: auto; /* required; upholds line height and height specification (where padding and line height determine overall height) */
  letter-spacing: -0.1px;
  line-height: 18px;
  padding: 6px 8px;
`;

/* Basic form field specification. */
/* Form label styles with targeted shared styles for form field elements i.e. form group, input, textarea and adornments (danger icon). */
export const StyledFormLabel = styled(Label)`
  /* Form label */
  &.${Classes.LABEL} {
    margin-bottom: 0;
  }

  /* Form group */
  .${Classes.FORM_GROUP} {
    margin-bottom: 0;
  }

  /* Input (shared styles for input and textarea) */
  &.${Classes.LABEL} .${Classes.INPUT} {
    margin-top: 8px; /* required; overrides BP label input margin top */
  }

  /* Input (shared styles for input and textarea) */
  .${Classes.INPUT} {
    ${formField};

    &:focus {
      border: 1px solid ${PRIMARY_BLUE} !important; /* required; overrides global.scss input border specification with important style declaration */
    }
  }

  /* Text area */
  textarea.${Classes.INPUT} {
    height: 32px; /* required; used by TextArea component to set the height of the component on first mount */
    min-height: 32px;
    min-width: 100%;
  }

  /* Danger input */
  textarea.${Classes.INPUT}.${Classes.INTENT_DANGER},
    .${Classes.INTENT_DANGER}
    .${Classes.INPUT} {
    border: 1px solid ${RED.C} !important; /* required; overrides global.scss input border specification with important style declaration */

    &:focus {
      box-shadow: none;
    }
  }

  /* Danger icon */
  .${Classes.ICON}.${Classes.INTENT_DANGER} {
    color: ${RED.C};
  }

  /* Danger icon positioner for input field */
  .${Classes.INPUT_GROUP} .${Classes.INPUT_ACTION}:last-child {
    right: 8px;
    top: 8px;
  }

  /* Adornment text */
  .${Classes.INPUT_GROUP} .${Classes.INPUT_LEFT_CONTAINER} {
    top: 0;
  }

  /* Helper text */
  .${Classes.FORM_HELPER_TEXT} {
    line-height: 15px;
    margin-top: 8px;
  }

  .${Classes.INTENT_DANGER} .${Classes.FORM_HELPER_TEXT} {
    color: ${RED.C};
  }
`;

/* Basic form field select button specification. */
/* Form label styles for select "button" (button to resemble basic form field specification). */
export const SelectFormLabel = styled.div`
  .${Classes.BUTTON} {
    ${formField};
    border: 1px solid ${LIGHT_GRAY.A}; /* mimics basic form field specification */
    display: flex;
    justify-content: space-between;
    margin-top: 8px; /* mimics basic form field specification */

    &:focus {
      outline: none;
    }
  }
`;

export const FormLabelText = styled.span`
  color: ${PT_TEXT_COLOR};
  display: block; /* required for form select "button" */
  letter-spacing: -0.1px;
  line-height: 18px;

  /* Optional label */
  > i {
    color: ${GRAY.A};
    padding-left: 4px;
  }
`;
