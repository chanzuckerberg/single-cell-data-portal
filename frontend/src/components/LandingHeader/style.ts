import { Classes, Colors } from "@blueprintjs/core";
import styled from "@emotion/styled";
import { Chip } from "czifui";
import {
  GRAY,
  PRIMARY_BLUE,
  PT_GRID_SIZE_PX,
  PT_TEXT_COLOR,
} from "../common/theme";

export const HEADER_HEIGHT_PX = 48;

export const Wrapper = styled.div`
  background-color: ${PT_TEXT_COLOR};
  height: ${HEADER_HEIGHT_PX}px;
  position: fixed;
  top: 0;
  width: 100%;
  z-index: 2;

  @media (max-width: 768px) {
    height: 100%;
    position: relative;
    display: flex;
  }
`;

export const MainWrapper = styled.div`
  align-items: center;
  display: flex;
  height: inherit; /* Take up full height of parent. */
  justify-content: space-between;
  padding: 0 16px;

  @media (max-width: 768px) {
    flex-direction: column;
    justify-content: flex-start;
    padding-top: ${HEADER_HEIGHT_PX}px;
  }
`;

export const MobileHomeLink = styled.span`
  display: none;
  @media (max-width: 768px) {
    display: block;
  }
`;

export const DesktopHomeLink = styled.span`
  @media (max-width: 768px) {
    display: none;
  }
`;

export const MobileMenuButton = styled.div`
  display: none;

  @media (max-width: 768px) {
    display: block;
    color: white;
    cursor: pointer;
    position: relative;
    z-index: 3;
    overflow: hidden;
  }
`;

export const MobileMenuButtonBar = styled.div`
  display: none;

  @media (max-width: 768px) {
    display: block;
    width: 30px;
    height: 2px;
    background-color: #fff;
    margin: 6px 0;
    transition: 0.4s;

    &.open {
      &:first-of-type {
        transform: translateY(8px) rotate(45deg);
      }

      &:nth-of-type(2) {
        transform: translateX(60px);
      }

      &:nth-of-type(3) {
        transform: translateY(-8px) rotate(-45deg);
      }
    }
  }
`;

export const MobileNavTray = styled.div`
  @media (max-width: 768px) {
    position: absolute;
    transform: translateX(100vw);
    right: 0;
    top: 0;
    width: 80%;
    height: 100vh;
    transition: transform 0.35s ease-out;

    &.active {
      transform: translateX(0);
    }
  }
`;

export const MobileNavWrapper = styled.div`
  position: relative;
  background-color: ${PT_TEXT_COLOR};
  width: 100%;
  height: ${HEADER_HEIGHT_PX}px;

  @media (max-width: 768px) {
    position: fixed;
    z-index: 3;
    top: 0;
    display: flex;
    justify-content: space-between;
    align-items: center;
    padding: 0 16px;
  }
`;

export const Left = styled.span`
  align-items: center;
  display: flex;
  gap: 32px;

  a {
    display: flex; /* Ensures the anchor wrapping the logo has correct line height. */
  }

  @media (max-width: 768px) {
    flex-direction: column;
    align-items: flex-start;
  }
`;

export const Right = styled.span`
  align-items: center;
  display: flex;
  gap: 24px;

  @media (max-width: 768px) {
    flex-direction: column;
    align-items: flex-start;
  }
`;

export const Nav = styled.span`
  display: flex;
  gap: 16px;

  @media (max-width: 768px) {
    flex-direction: column;
    align-items: flex-start;
    padding-top: 16px;
    position: relative;
    left: -14px;
    margin-bottom: 26px;
  }
`;

function button(): string {
  return `
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
}

function iconButton(): string {
  return `
    ${button}

    .${Classes.ICON} {
      color: inherit; /* Overrides BP button icon color rule by inheriting color from parent. */
    }
  `;
}

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

export const HiringLink = styled.a`
  background-color: ${PRIMARY_BLUE};
  padding: 7px 14px;
  border-radius: 4px;
  color: #fff;
  font-weight: 600;
  transition: 0.3s;
  &:hover {
    color: #fff;
    text-decoration: none !important;
    background-color: #0056c6;
  }
`;
