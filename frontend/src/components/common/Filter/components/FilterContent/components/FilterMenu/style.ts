import styled from "@emotion/styled";
import { CommonThemeProps } from "@czi-sds/components";
import {
  listCss,
  listItemButtonCss,
  listItemCss,
  listItemDividerCss,
  listItemIconCss,
  listItemTextCss,
  scrollbar,
} from "src/components/common/Filter/components/FilterContent/components/common/style";
import { spacesXs } from "src/common/theme";

export const MAX_DISPLAYABLE_MENU_ITEMS = 9;
const DIVIDER_HEIGHT_PX = 9;
const MENU_ITEM_HEIGHT_PX = 32;
const MAX_MENU_HEIGHT_PX =
  (MAX_DISPLAYABLE_MENU_ITEMS + 0.5) *
  MENU_ITEM_HEIGHT_PX; /* Undivided menu - max height is 9.5 menu items. */
const MAX_DIVIDED_MENU_HEIGHT_PX =
  MAX_MENU_HEIGHT_PX +
  DIVIDER_HEIGHT_PX; /* Divided menu - max undivided menu height, plus divider height. */

interface Props extends CommonThemeProps {
  isMenuDivided: boolean;
  menuWidth: number;
  scrollable: boolean;
}

export const FilterMenu = styled.span<Props>`
  min-width: ${(props) =>
    `${props.menuWidth}px`}; /* maintains menu width when filtering menu items */
  padding: ${spacesXs}px;

  .MuiList-root {
    ${listCss}
    ${scrollbar}
    max-height: ${({ isMenuDivided }) =>
      isMenuDivided
        ? `${MAX_DIVIDED_MENU_HEIGHT_PX}px`
        : `${MAX_MENU_HEIGHT_PX}px`};
    overflow-y: auto;
    padding-right: ${(props) =>
      props.scrollable ? `${spacesXs(props)}px` : undefined};

    li {
      ${listItemCss}
      ${listItemButtonCss}
      display: flex;
      cursor: pointer;

      /* "No Matches" menu item */

      &.MuiListItem-root {
        cursor: default;

        &:hover {
          background-color: transparent;
        }
      }

      .MuiListItemIcon-root {
        ${listItemIconCss}
      }

      .MuiListItemText-root {
        ${listItemTextCss}
      }
    }

    .MuiDivider-root {
      ${listItemDividerCss}
    }
  }
`;
