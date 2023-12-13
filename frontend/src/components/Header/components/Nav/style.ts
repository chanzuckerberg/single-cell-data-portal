import styled from "@emotion/styled";
import { Classes } from "@blueprintjs/core";
import { css } from "@emotion/react";
import {
  fontWeightBold,
  fontWeightSemibold,
  gray300,
  gray500,
  grayWhite,
  spacesL,
} from "src/common/theme";
import { CommonThemeProps, fontCapsXxxs } from "@czi-sds/components";

export const Nav = styled.span`
  display: flex;
  gap: ${spacesL}px;
`;

export const button = (props: CommonThemeProps) => css`
  display: inline-block; /* Wrapper to mimic line height of children. */

  .${Classes.BUTTON}.${Classes.MINIMAL} {
    background: none;
    border-radius: 0;
    color: ${gray300(props)};
    font-size: 13px;
    font-weight: ${fontWeightSemibold(props)};
    height: 22px;
    letter-spacing: -0.1px;
    line-height: 18px;
    min-height: 22px;
    padding: 0;

    &.${Classes.ACTIVE}, &:hover {
      color: ${grayWhite()};
    }

    &:focus {
      outline: none;
    }
  }
`;

export const LinkWrapper = styled.span`
  ${button}
  .${Classes.BUTTON}.${Classes.MINIMAL}.${Classes.ACTIVE} {
    box-shadow: inset 0 -2px 0 ${grayWhite} !important; /* Overrides specificity of BP button active box shadow rule. */
  }
  align-items: center;
  display: flex;
`;

export const NavSection = styled.span`
  display: flex;
  flex-direction: column;
  align-items: baseline;
`;

export const NavSectionTitle = styled.span`
  ${fontCapsXxxs}
  color: ${gray500};
  font-weight: ${fontWeightBold};
`;

export const NavItemContainer = styled.span`
  display: flex;
  gap: ${spacesL}px;
`;
