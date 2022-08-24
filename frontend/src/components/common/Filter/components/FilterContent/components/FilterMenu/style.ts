import { Classes, MenuItem } from "@blueprintjs/core";
import { Divider } from "@material-ui/core";
import { GRAY, PT_TEXT_COLOR } from "src/components/common/theme";
import styled from "styled-components";
import { scrollbar } from "../../../../common/style";

export const MAX_DISPLAYABLE_MENU_ITEMS = 9;
const DIVIDER_HEIGHT_PX = 9;
const MENU_ITEM_HEIGHT_PX = 32;
const MAX_MENU_HEIGHT_PX =
  (MAX_DISPLAYABLE_MENU_ITEMS + 0.5) *
  MENU_ITEM_HEIGHT_PX; /* Undivided menu - max height is 9.5 menu items. */
const MAX_DIVIDED_MENU_HEIGHT_PX =
  MAX_MENU_HEIGHT_PX +
  DIVIDER_HEIGHT_PX; /* Divided menu - max undivided menu height, plus divider height. */

interface Props {
  menuWidth: number;
}

interface MenuListProps {
  isMenuDivided?: boolean;
}

export const MenuWrapper = styled.span<Props>`
  .${Classes.MENU} {
    min-width: ${(props) =>
      `${props.menuWidth}px`}; /* overrides BP menu min-width specification; maintains menu width when filtering menu items */
    padding: 6px;

    li {
      margin: 0; /* overrides margin from layout.css */

      &:focus {
        outline: none;
      }
    }
  }
`;

export const MenuItemsWrapper = styled.div<MenuListProps>`
  max-height: ${({ isMenuDivided = false }) =>
    isMenuDivided
      ? `${MAX_DIVIDED_MENU_HEIGHT_PX}px`
      : `${MAX_MENU_HEIGHT_PX}px`};
  overflow-y: auto;
  padding-right: 6px;
  ${scrollbar};
`;

export const MenuDivider = styled(Divider)`
  background-color: ${PT_TEXT_COLOR};
  opacity: 0.15;
  margin: 4px 0;
`;

export const NoMatches = styled(MenuItem)`
  && {
    color: ${GRAY.A};
    letter-spacing: -0.1px;
    line-height: 18px;
    padding: 7px 8px;

    &:hover {
      background-color: transparent;
      color: ${GRAY.A};
      cursor: default;
    }
  }
`;
