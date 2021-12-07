import { Classes, MenuItem } from "@blueprintjs/core";
import { GRAY, LIGHT_GRAY, PT_TEXT_COLOR } from "src/components/common/theme";
import styled from "styled-components";

export const MAX_DISPLAYABLE_MENU_ITEMS = 9;

interface Props {
  menuWidth: number;
}

interface MenuItemProps {
  isSelected: boolean;
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

export const InputGroupWrapper = styled.div`
  margin-bottom: 4px;

  .${Classes.INPUT_GROUP} {
    .${Classes.ICON} {
      color: ${GRAY.A};
      margin: 0;
      padding: 4px;
      top: 50%;
      transform: translateY(-50%);
    }

    .${Classes.INPUT} {
      border: 1px solid ${LIGHT_GRAY.A} !important; /* required; overrides BP input border with important style declaration */
      border-radius: 3px;
      color: ${PT_TEXT_COLOR};
      letter-spacing: -0.1px;
      line-height: 18px;
      height: 32px;

      &::placeholder {
        opacity: 0.6;
      }
    }
  }
`;

export const MenuItemsWrapper = styled.div`
  max-height: calc((${MAX_DISPLAYABLE_MENU_ITEMS} + 0.5) * 32px);
  overflow-y: auto;
  padding-right: 6px;

  &::-webkit-scrollbar {
    width: 4px;
  }

  &::-webkit-scrollbar-thumb {
    background-clip: content-box;
    background-color: ${LIGHT_GRAY.A};
    border-radius: 4px;
  }
`;

export const NoMatches = styled(MenuItem)`
  color: ${GRAY.A};
  letter-spacing: -0.1px;
  line-height: 18px;

  &:hover {
    background-color: transparent;
    color: ${GRAY.A};
    cursor: default;
  }
`;

export const MenuItemWrapper = styled.span<MenuItemProps>`
  .${Classes.MENU_ITEM} {
    color: ${PT_TEXT_COLOR};
    letter-spacing: -0.1px;
    line-height: 18px;
    padding: 7px 8px;

    .${Classes.FILL} {
      font-weight: ${(props: MenuItemProps) =>
        props.isSelected ? 500 : undefined};
      margin-right: 8px;
    }

    &:focus {
      outline: none;
    }

    .${Classes.ICON} {
      align-items: center;
      color: #0073ff;
      display: flex;
      height: 18px;
      justify-content: center;
      margin: 0 8px 0 0;

      > svg {
        height: 14px;
        width: 14px;
      }
    }
  }
`;
