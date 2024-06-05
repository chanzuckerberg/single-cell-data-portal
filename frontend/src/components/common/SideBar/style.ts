import styled from "@emotion/styled";
import { LIGHT_GRAY } from "src/components/common/theme";
import { HEADER_HEIGHT_PX } from "src/components/Header/style";
import { CommonThemeProps, fontHeaderM } from "@czi-sds/components";
import {
  grey300,
  grey500,
  spacesL,
  spacesS,
  spacesXl,
  spacesXxxs,
} from "src/common/theme";
import { css, SerializedStyles } from "@emotion/react";
import { Button } from "src/components/common/Button";

export enum Position {
  LEFT = "left",
  RIGHT = "right",
}

interface Props extends CommonThemeProps {
  sideBarWidth: number;
  position: (typeof Position)[keyof typeof Position];
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
  background-color: white;
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

export const ToggleButtonText = styled.span`
  ${fontHeaderM}
  color: #000000;
`;

export const ToggleButton = styled(Button, {
  shouldForwardProp: (prop) => prop !== "isExpanded",
})<PositionerProps>`
  .MuiButton-endIcon {
    margin: 0;
  }

  &:hover {
    .MuiButton-endIcon {
      .MuiSvgIcon-root {
        color: ${grey500};
      }
    }
  }

  ${(props) =>
    props.isExpanded &&
    css`
      justify-content: space-between;
      margin-bottom: ${spacesS(props)}px;
      padding: 0;
    `}
  ${(props) =>
    !props.isExpanded &&
    css`
      flex-direction: column-reverse;
      gap: ${spacesS(props)}px;
      min-width: 0;
      padding: ${spacesXl(props)}px ${spacesS(props)}px;

      ${ToggleButtonText} {
        transform: scale(-1);
        writing-mode: vertical-rl;
      }
    `}
`;

export function sideBarPositionerPadding(
  props: PositionerProps
): SerializedStyles | undefined {
  if (props.isExpanded) {
    return css`
      padding: ${spacesL(props)}px ${spacesS(props)}px;

      ${ToggleButton} {
        margin-bottom: ${spacesXxxs(props)}px;
        padding: 0 ${spacesS(props)}px;
      }
    `;
  }
}
