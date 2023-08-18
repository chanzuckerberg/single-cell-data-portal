import styled from "@emotion/styled";
import { LIGHT_GRAY } from "src/components/common/theme";
import { HEADER_HEIGHT_PX } from "src/components/Header/style";
import { CommonThemeProps, fontBodyS } from "@czi-sds/components";
import {
  grey300,
  grey500,
  fontWeightSemibold,
  spacesL,
  spacesS,
  spacesXl,
} from "src/common/theme";

export enum Position {
  LEFT = "left",
  RIGHT = "right",
}

interface Props extends CommonThemeProps {
  sideBarWidth: number;
  position: typeof Position[keyof typeof Position];
}

interface PositionerProps extends CommonThemeProps {
  isExpanded: boolean;
}

const generateBoxShadow = (props: Props) =>
  `inset ${props.position === Position.LEFT ? -0.5 : 0.5}px 0 ${grey300(
    props
  )}`;

export const SideBar = styled.div<Props>`
  box-sizing: border-box;
  box-shadow: ${generateBoxShadow};
  width: ${(props) => `${props.sideBarWidth}px`};
  grid-area: ${(props) =>
    props.position === Position.LEFT ? "leftsidebar" : "rightsidebar"};
`;

export const SideBarPositioner = styled.div<PositionerProps>`
  max-height: calc(
    100vh - ${HEADER_HEIGHT_PX}px
  ); /* required for sidebar scrolling where header height is 48px */
  overflow-y: ${(props) =>
    props.isExpanded ? "overlay" : undefined}; /* overlay is deprecated */
  padding: ${(props) =>
    props.isExpanded ? `${spacesXl(props)}px ${spacesL(props)}px` : undefined};
  position: sticky;
  top: ${HEADER_HEIGHT_PX}px;
  width: inherit; /* inherits sidebar container width specification */

  &::-webkit-scrollbar {
    width: 16px; /* accommodates track size 4px and the -webkit-scrollbar-thumb 6px transparent border */
  }

  &::-webkit-scrollbar-thumb {
    background-clip: content-box;
    background-color: ${LIGHT_GRAY.A};
    border: 6px solid transparent;
    border-radius: 12px;
  }
`;

export const SideBarToggleButtonWrapper = styled.span`
  display: block;

  .MuiButton-root {
    width: 100%;

    &:hover {
      background: none;
    }

    .MuiButton-endIcon {
      align-items: center;
      display: flex;
      height: 20px;
      justify-content: center;
      margin: 0;
      width: 20px;

      .MuiSvgIcon-root {
        color: ${grey500};
        height: 14px;
        width: 14px;
      }
    }
  }
`;

export const ToggleButtonText = styled.span`
  ${fontBodyS}
  color: #000000;
  font-weight: ${fontWeightSemibold};
  letter-spacing: -0.006em;
`;

export const SideBarClosedButtonWrapper = styled(SideBarToggleButtonWrapper)`
  .MuiButton-root {
    flex-direction: column-reverse;
    gap: ${spacesS}px;
    min-width: 0;
    padding: ${spacesXl}px ${spacesS}px;

    ${ToggleButtonText} {
      transform: scale(-1);
      writing-mode: vertical-rl;
    }
  }
`;

export const SideBarOpenButtonWrapper = styled(SideBarToggleButtonWrapper)`
  .MuiButton-root {
    justify-content: space-between;
    margin-bottom: ${spacesS}px;
    padding: 0;
  }
`;
