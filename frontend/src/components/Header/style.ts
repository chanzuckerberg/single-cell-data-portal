import { Classes, Colors } from "@blueprintjs/core";
import { css } from "@emotion/react";
import styled from "@emotion/styled";
import { Chip, getColors, InputDropdown } from "@czi-sds/components";
import { PT_GRID_SIZE_PX, PT_TEXT_COLOR } from "../common/theme";
import { button } from "src/components/Header/components/Nav/style";
import { spacesL, spacesXl } from "src/common/theme";

export const HEADER_HEIGHT_PX = 48;

export const Wrapper = styled.div`
  background-color: ${PT_TEXT_COLOR};
  height: ${HEADER_HEIGHT_PX}px;
  position: fixed;
  top: 0;
  width: 100%;
  z-index: 1;
`;

export const MainWrapper = styled.div`
  align-items: center;
  display: flex;
  gap: ${spacesXl}px;
  height: inherit; /* Take up full height of parent. */
  justify-content: space-between;
  padding: 0 16px;
`;

export const Left = styled.span`
  align-items: center;
  display: flex;
  gap: 32px;

  a {
    display: flex; /* Ensures the anchor wrapping the logo has correct line height. */
  }
`;

export const Right = styled.span`
  align-items: center;
  display: flex;
  gap: ${spacesL}px;
`;

const iconButton = css`
  ${button}
  .${Classes.ICON} {
    color: inherit; /* Overrides BP button icon color rule by inheriting color from parent. */
  }
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

export const StyledInputDropdown = styled(InputDropdown)`
  /* Overriding colors after SDS v14.5.0 upgrade */
  ${(props) => {
    const colors = getColors(props);

    return `
      span {
        color: ${colors?.gray["500"]} !important;
      }

      :hover {
        background: none;
        span {
          color: ${colors?.gray["600"]} !important;
        }
      }
    `;
  }}
`;
