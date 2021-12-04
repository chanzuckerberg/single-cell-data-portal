import { Classes } from "@blueprintjs/core";
import { GRAY } from "src/components/common/theme";
import styled from "styled-components";

export const CategoryButton = styled.span`
  display: block;

  .${Classes.BUTTON} {
    color: ${GRAY.A};
    display: flex;
    font-weight: 500;
    justify-content: flex-start;
    letter-spacing: -0.1px;
    line-height: 18px;
    min-height: 20px; /* overrides specificity of bp3 button min height rule */
    height: 20px; /* overrides specificity of bp3 button height rule */
    padding: 0;
    width: 100%;

    &:hover {
      color: ${GRAY.A};
      background: none;
    }

    &:focus {
      outline: none;
    }

    .${Classes.BUTTON_TEXT} {
      height: inherit; /* required; inheriting button height accurately center aligns button text. */
    }

    &.${Classes.MINIMAL}:disabled {
      color: ${GRAY.E};
      cursor: auto;
    }
  }
`;
