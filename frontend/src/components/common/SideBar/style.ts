import { Classes } from "@blueprintjs/core";
import { LIGHT_GRAY, PT_TEXT_COLOR } from "src/components/common/theme";
import styled from "styled-components";

interface Props {
  sideBarWidth: number;
}

interface PositionerProps {
  isExpanded: boolean;
}

export const SideBar = styled.div<Props>`
  box-sizing: border-box;
  box-shadow: inset -1px 0px 0px rgba(16, 22, 26, 0.15);
  width: ${(props: Props) => `${props.sideBarWidth}px`};
`;

export const SideBarPositioner = styled.div<PositionerProps>`
  max-height: calc(
    100vh - 48px
  ); /* required for sidebar scrolling where header height is 48px */
  overflow-y: ${(props: PositionerProps) =>
    props.isExpanded ? "overlay" : undefined};
  padding: ${(props: PositionerProps) =>
    props.isExpanded ? "24px 16px" : undefined};
  position: sticky;
  top: 48px;
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

const SideBarToggleButtonWrapper = styled.span`
  .${Classes.BUTTON} {
    color: ${PT_TEXT_COLOR};
    display: flex;
    /* TODO(cc) font weight specification correct; font rendering resolves with https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/chanzuckerberg/single-cell-data-portal/1685 */
    font-weight: 500;
    letter-spacing: -0.1px;
    line-height: 18px;
    min-height: 20px; /* overrides specificity of bp3 button min height rule */

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
    height: 20px; /* overrides specificity of bp3 button height rule */
    justify-content: space-between;
    padding: 0;
    width: 100%;
  }
`;
