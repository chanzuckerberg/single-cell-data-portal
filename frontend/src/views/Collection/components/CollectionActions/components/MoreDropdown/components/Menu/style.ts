import styled from "@emotion/styled";
import {
  CommonThemeProps,
  IconNameToSmallSizes,
  Menu as SDSMenu,
  MenuItem as SDSMenuItem,
  MenuItemProps as SDSMenuItemProps,
} from "@czi-sds/components";
import { error400, gray400, shadowM, spacesM, spacesS } from "src/common/theme";
import { MENU_ITEM_COLOR } from "src/views/Collection/components/CollectionActions/components/MoreDropdown/components/Menu/types";
import { css, SerializedStyles } from "@emotion/react";

export interface StyledMenuItemProps<
  IconName extends keyof IconNameToSmallSizes,
> extends SDSMenuItemProps<IconName>,
    CommonThemeProps {
  color?: MENU_ITEM_COLOR;
}

export const StyledMenu = styled(SDSMenu)`
  &.MuiPopover-root {
    z-index: 10; /* positions menu beneath BP portals to allow proper layering with delete/publish collection alert interactions */
  }

  .MuiMenu-paper {
    box-shadow: ${shadowM};
    min-width: 199px; /* matches mocks */
  }
`;

export const StyledMenuItem = styled(SDSMenuItem)<
  StyledMenuItemProps<"TrashCan" | "Edit">
>`
  &.MuiMenuItem-root {
    ${customMenuItem}
    .primary-text {
      > span {
        display: flex;
        margin-right: ${spacesS}px;
      }

      svg {
        align-self: center;
      }
    }
  }
`;

function customMenuItem(
  props: StyledMenuItemProps<"TrashCan" | "Edit">
): SerializedStyles | undefined {
  // Custom "error" color for menu item (with "error" icon).
  if (props.color === MENU_ITEM_COLOR.ERROR) {
    return css`
      &:not(.Mui-disabled) {
        .primary-text {
          color: ${error400(props)};
        }
      }
    `;
  }
  // Custom "gray" color for menu item (with custom icon).
  if (props.color === MENU_ITEM_COLOR.GRAY) {
    return css`
      .MuiSvgIcon-root {
        color: ${gray400(props)};
        font-size: ${spacesM(props)}px;
      }

      &.Mui-disabled {
        .MuiSvgIcon-root {
          color: inherit;
        }
      }
    `;
  }
}
