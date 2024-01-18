import styled from "@emotion/styled";
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

export interface ButtonProps extends CommonThemeProps {
  isActive?: boolean;
}

export const button = (props: ButtonProps) => css`
  display: inline-block; /* Wrapper to mimic line height of children. */

  background: none;
  border-radius: 0;

  font-size: 13px;
  font-weight: ${fontWeightSemibold(props)};
  height: 22px;
  letter-spacing: -0.1px;
  line-height: 18px;
  min-height: 22px;
  padding: 0;

  a {
    color: ${props.isActive ? grayWhite() : gray300(props)};

    text-decoration: ${props.isActive ? "underline white solid 2px" : "none"};
    text-underline-offset: 4px;

    &:hover {
      color: ${grayWhite()};
    }

    &:focus {
      outline: none;
    }

    &:focus-visible {
      outline: 1px solid white;
      outline-offset: 3px;
    }
  }
`;

export const LinkWrapper = styled.span`
  ${button}

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
