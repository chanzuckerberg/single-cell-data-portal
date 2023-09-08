import styled from "@emotion/styled";
import { css } from "@emotion/react";
import { IconButton } from "@mui/material";

export const FullscreenButton = styled(IconButton)`
  visibility: hidden;
  transition: visibility 0.2s;
  z-index: 2;
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
          z-index: 9999;
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
