import styled from "@emotion/styled";
import { css } from "@emotion/react";
import { fontBodyXxxs, getColors } from "czifui";
import { IconButton } from "@mui/material";

export const StyledLegendText = styled.text`
  ${fontBodyXxxs}

  alignment-baseline: middle;

  ${(props) => {
    const colors = getColors(props);
    return `color: ${colors?.gray[500]}`;
  }}
`;

export const FullscreenButton = styled(IconButton)`
  visibility: hidden;
  transition: visibility 0.2s;
  z-index: 1;
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
