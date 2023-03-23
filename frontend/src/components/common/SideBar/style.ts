import { Classes } from "@blueprintjs/core";
import styled from "@emotion/styled";
import { LIGHT_GRAY, PT_TEXT_COLOR } from "src/components/common/theme";
import { HEADER_HEIGHT_PX } from "src/components/Header/style";

export enum Position {
  LEFT = "left",
  RIGHT = "right",
}

interface Props {
  sideBarWidth: number;
  position: typeof Position[keyof typeof Position];
}

interface PositionerProps {
  isExpanded: boolean;
}

export const SIDEBAR_BOX_SHADOW_COLOR = "rgba(16, 22, 26, 0.15)";

const generateBoxShadow = (props: Props) =>
  `inset ${
    props.position === Position.LEFT ? -1 : 1
  }px 0px 0px ${SIDEBAR_BOX_SHADOW_COLOR}`;

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
  ); /* required for sidebar scrolling  where header height is 48px */
  overflow-y: ${(props) =>
    props.isExpanded ? "overlay" : undefined}; /* overlay is deprecated */
  padding: ${(props) => (props.isExpanded ? "24px 16px" : undefined)};
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
  .${Classes.BUTTON} {
    color: ${PT_TEXT_COLOR};
    display: flex;
    font-weight: 500;
    letter-spacing: -0.1px;
    line-height: 18px;
    min-height: 20px; /* overrides specificity of bp4 button min height rule */

    &:hover {
      background: none;
    }

    &:focus {
      outline: none;
    }

    .${Classes.ICON} > svg {
      height: 20px;
      width: 20px;
    }
  }
`;

export const SideBarClosedButtonWrapper = styled(SideBarToggleButtonWrapper)`
  .${Classes.BUTTON} {
    flex-direction: column-reverse;
    gap: 12px;
    height: 100%;
    justify-content: flex-end;
    padding: 24px 8px;
    width: 100%; /* button click target consumes full width of closed sidebar */

    .${Classes.BUTTON_TEXT} {
      margin: 0;
      transform: scale(-1);
      writing-mode: vertical-rl;
    }
  }
`;

export const SideBarOpenButtonWrapper = styled(SideBarToggleButtonWrapper)`
  display: block;
  margin-bottom: 12px;

  .${Classes.BUTTON} {
    height: 20px; /* overrides specificity of bp4 button height rule */
    justify-content: space-between;
    padding: 0;
    width: 100%;
  }
`;
