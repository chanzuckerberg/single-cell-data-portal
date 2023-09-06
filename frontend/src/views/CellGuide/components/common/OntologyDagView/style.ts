import styled from "@emotion/styled";
import { css } from "@emotion/react";
import { IconButton } from "@mui/material";
import { HEADER_HEIGHT_PX } from "src/components/LandingHeader/style";

interface FullscreenButtonProps {
  isFullScreen: boolean;
}
export const FullscreenButton = styled(IconButton)<FullscreenButtonProps>`
  visibility: hidden;
  transition: visibility 0.2s;
  z-index: 99;
  margin-top: ${(props) => (props.isFullScreen ? HEADER_HEIGHT_PX : 0)}px;
`;

interface HoverContainerProps {
  isFullScreen: boolean;
  height: number;
  width: number;
}

export const HoverContainer = styled.div<HoverContainerProps>`
  position: relative;
  ${({ isFullScreen, height, width }) =>
    isFullScreen
      ? css`
          position: fixed;
          top: 0;
          left: 0;
          width: 100%;
          height: 100%;
          overflow: auto;
          z-index: 1;
        `
      : css`
          height: ${height}px;
          width: ${width}px;
        `}

  &:hover {
    ${FullscreenButton} {
      visibility: visible;
    }
  }
`;

export const TooltipInPortalStyle = css`
  .visx-tooltip {
    z-index: 1;
  }
`;

interface StyledSVGProps {
  isDragging: boolean;
}
export const StyledSVG = styled.svg<StyledSVGProps>`
  cursor: ${({ isDragging }) => (isDragging ? "grabbing" : "grab")};
  touch-action: none;
  position: absolute;
  top: 0;
  left: 0;
`;
