import { Classes, Colors } from "@blueprintjs/core";
import styled from "@emotion/styled";
import { css } from "@emotion/react";
import { GRAY, PT_TEXT_COLOR } from "./theme";
import { HEADER_HEIGHT_PX } from "../../globals";
import { Chip } from "czifui";

export const Wrapper = styled.div`
  background-color: ${PT_TEXT_COLOR};
  height: ${HEADER_HEIGHT_PX}px;
  position: fixed;
  top: 0;
  width: 100%;
  z-index: 1;
`;

export const MainWrapper = styled.div`
  align-items: center;
  display: flex;
  height: inherit; /* Take up full height of parent. */
  justify-content: space-between;
  padding: 0 16px;
`;

export const Left = styled.span`
  align-items: center;
  display: flex;
  gap: 32px;

  a {
    display: flex; /* Ensures the anchor wrapping the logo has correct line height. */
  }
`;

export const Right = styled.span`
  align-items: center;
  display: flex;
  gap: 24px;
`;

export const Nav = styled.span`
  display: flex;
  gap: 16px;
`;

const button = css`
  display: inline-block; /* Wrapper to mimic line height of children. */

  .${Classes.BUTTON}.${Classes.MINIMAL} {
    background: none;
    border-radius: 0;
    color: ${GRAY.D};
    height: 22px;
    letter-spacing: -0.1px;
    line-height: 18px;
    min-height: 22px;
    padding: 0;

    > span {
      font-size: 13px;
      font-weight: 500;
    }

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

  display: flex;
  align-items: center;
`;

export const BetaChip = styled(Chip)`
  background: #7a41ce;
  color: white;
  margin-left: 4px;
  height: 16px !important;
`;
