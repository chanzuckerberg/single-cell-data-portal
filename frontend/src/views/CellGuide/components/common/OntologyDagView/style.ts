import styled from "@emotion/styled";
import { css } from "@emotion/react";
import { IconButton } from "@mui/material";
import { ButtonIcon, TagFilter } from "@czi-sds/components";
import { primary400, spacesL } from "src/common/theme";
import { HEADER_HEIGHT_PX } from "src/components/Header/style";

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

export const RightAligned = styled.div`
  width: 100%;
  display: flex;
  flex-direction: row;
  justify-content: flex-end;
  padding-right: ${spacesL}px;
`;

export const StyledButtonIcon = styled(ButtonIcon)`
  z-index: 1;
`;

export const HoverContainer = styled.div<HoverContainerProps>`
  position: relative;
  ${({ isFullScreen, height, width }) =>
    isFullScreen
      ? css`
          position: fixed;
          top: ${HEADER_HEIGHT_PX}px;
          left: 0;
          width: ${width}px;
          height: ${height}px;
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
    z-index: 99999;
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

export const StyledTagFilter = styled(TagFilter)`
  z-index: 1;
  .MuiChip-label {
    color: ${primary400};
  }
  .MuiSvgIcon-root {
    fill: ${primary400};
  }
  &:hover {
    background-color: #e0f0ff;
    .MuiChip-label {
      color: ${primary400};
    }
    .MuiSvgIcon-root {
      fill: ${primary400};
    }
  }
`;
