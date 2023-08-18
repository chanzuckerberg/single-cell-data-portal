import styled from "@emotion/styled";
import { Classes, Colors } from "@blueprintjs/core";
import { GRAY } from "src/components/common/theme";
import { css } from "@emotion/react";
import { spacesL } from "src/common/theme";

export const Nav = styled.span`
  display: flex;
  gap: ${spacesL}px;
`;

export const button = css`
  display: inline-block; /* Wrapper to mimic line height of children. */

  .${Classes.BUTTON}.${Classes.MINIMAL} {
    background: none;
    border-radius: 0;
    color: ${GRAY.D};
    font-size: 13px;
    font-weight: 500;
    height: 22px;
    letter-spacing: -0.1px;
    line-height: 18px;
    min-height: 22px;
    padding: 0;

    &.${Classes.ACTIVE}, &:hover {
      color: ${Colors.WHITE};
    }

    &:focus {
      outline: none;
    }
  }
`;

export const LinkWrapper = styled.span`
  ${button}
  .${Classes.BUTTON}.${Classes.MINIMAL}.${Classes.ACTIVE} {
    box-shadow: inset 0 -2px 0 ${Colors.WHITE} !important; /* Overrides specificity of BP button active box shadow rule. */
  }

  align-items: center;
  display: flex;
`;
