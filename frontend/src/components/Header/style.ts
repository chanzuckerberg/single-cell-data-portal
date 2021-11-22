import { Classes, Colors } from "@blueprintjs/core";
import styled, { css } from "styled-components";
import { layout } from "../common/layout";
import { GRAY, PT_TEXT_COLOR } from "../common/theme";

export const Wrapper = styled.div`
  background-color: ${PT_TEXT_COLOR};
  height: 48px;
  margin-bottom: 24px;
  position: sticky;
  top: 0;
  width: 100%;
  z-index: 1;
`;

export const MainWrapper = styled.div`
  ${layout};
  align-items: center;
  display: flex;
  height: inherit; /* Take up full height of parent. */
  justify-content: space-between;
  margin: 0 auto;
  padding: 0 16px;
`;

export const Left = styled.span`
  a {
    display: flex; /* Ensures the anchor wrapping the logo has correct line height. */
  }
`;

export const Right = styled.span`
  display: flex;
  gap: 24px;
`;

const button = css`
  display: inline-block; /* Wrapper to mimic line height of children. */

  .${Classes.BUTTON}.${Classes.MINIMAL} {
    background: none;
    border-radius: 0;
    color: ${GRAY.D};
    font-size: 12px;
    font-weight: 600;
    height: 23px;
    line-height: 15px;
    min-height: 23px;
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
  ${button};

  .${Classes.BUTTON}.${Classes.MINIMAL}.${Classes.ACTIVE} {
    box-shadow: inset 0 -2px 0 ${Colors.WHITE} !important; /* Overrides specificity of BP button active box shadow rule. */
  }
`;

export const LearnButtonWrapper = styled.span`
  ${button};
`;

export const AuthButtonWrapper = styled.span`
  ${button};

  .${Classes.BUTTON}.${Classes.MINIMAL} {
    font-weight: 400; /* Overrides locally defined button font weight rule. */
  }

  .${Classes.BUTTON}.${Classes.MINIMAL}, .${Classes.ICON} {
    color: ${Colors.WHITE}; /* Overrides locally defined button and BP button icon color rule. */
  }
`;
