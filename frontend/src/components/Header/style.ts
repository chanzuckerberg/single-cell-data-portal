import { Classes, Colors } from "@blueprintjs/core";
import { Chip } from "czifui";
import styled, { css } from "styled-components";
import { GRAY, PT_GRID_SIZE_PX, PT_TEXT_COLOR } from "../common/theme";

export const HEADER_HEIGHT_PX = 48;

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

const iconButton = css`
  ${button}

  .${Classes.ICON} {
    color: inherit; /* Overrides BP button icon color rule by inheriting color from parent. */
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

export const LearnButtonWrapper = styled.span`
  ${iconButton}
`;

export const AuthButtonWrapper = styled.span`
  ${iconButton}

  .${Classes.BUTTON}.${Classes.MINIMAL} {
    color: ${Colors.WHITE}; /* Overrides locally defined button color rule. */
    font-weight: 400; /* Overrides locally defined button font weight rule. */
  }
`;

export const BetaChip = styled(Chip)`
  background: #7a41ce;
  color: white;
  margin-left: ${PT_GRID_SIZE_PX / 2}px;
  height: ${PT_GRID_SIZE_PX * 2}px !important;
`;
