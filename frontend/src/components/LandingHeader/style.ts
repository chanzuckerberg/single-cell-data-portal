import styled from "@emotion/styled";
import { PT_TEXT_COLOR } from "../common/theme";
import Nav from "src/components/Header/components/Nav";
import { NavDivider } from "src/components/Header/components/Nav/components/NavDivider/style";
import { spacesL, spacesXl } from "src/common/theme";
import { fontHeaderM } from "@czi-sds/components";

export const HEADER_HEIGHT_PX = 48;
export const MAX_WIDTH_BREAKPOINT_PX = 768;

export const Wrapper = styled.div`
  background-color: ${PT_TEXT_COLOR};
  height: ${HEADER_HEIGHT_PX}px;
  position: fixed;
  top: 0;
  width: 100%;
  z-index: 2;

  @media (max-width: ${MAX_WIDTH_BREAKPOINT_PX}px) {
    height: 100%;
    position: relative;
    display: flex;
  }
`;

export const MainWrapper = styled.div`
  align-items: center;
  display: flex;
  gap: ${spacesXl}px;
  height: inherit; /* Take up full height of parent. */
  justify-content: space-between;
  padding: 0 16px;

  @media (max-width: ${MAX_WIDTH_BREAKPOINT_PX}px) {
    align-items: flex-start;
    flex-direction: column;
    justify-content: flex-start;
  }
`;

export const MobileHomeLink = styled.span`
  display: none;
  @media (max-width: ${MAX_WIDTH_BREAKPOINT_PX}px) {
    display: block;
  }
`;

export const DesktopHomeLink = styled.span`
  @media (max-width: ${MAX_WIDTH_BREAKPOINT_PX}px) {
    display: none;
  }
`;

export const MobileMenuButton = styled.div`
  display: none;

  @media (max-width: ${MAX_WIDTH_BREAKPOINT_PX}px) {
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

  @media (max-width: ${MAX_WIDTH_BREAKPOINT_PX}px) {
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
  @media (max-width: ${MAX_WIDTH_BREAKPOINT_PX}px) {
    position: absolute;
    transform: translateX(100vw);
    right: 0;
    top: ${HEADER_HEIGHT_PX}px;
    width: 100%;
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

  @media (max-width: ${MAX_WIDTH_BREAKPOINT_PX}px) {
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

  @media (max-width: ${MAX_WIDTH_BREAKPOINT_PX}px) {
    flex-direction: column;
    align-items: flex-start;
  }
`;

export const Right = styled.span`
  align-items: center;
  display: flex;
  gap: ${spacesL}px;

  @media (max-width: ${MAX_WIDTH_BREAKPOINT_PX}px) {
    flex-direction: column;
    align-items: flex-start;
  }
`;

export const Navigation = styled(Nav)`
  @media (max-width: ${MAX_WIDTH_BREAKPOINT_PX}px) {
    flex-direction: column;
    align-items: flex-start;
    padding-top: 16px;

    ${NavDivider} {
      height: 1px;
      width: 100%;
    }
  }
`;

export const HeaderTitle = styled.div`
  ${fontHeaderM}
  color: white;
  cursor: pointer;
`;
