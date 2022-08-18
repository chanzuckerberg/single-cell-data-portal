import { Classes } from "@blueprintjs/core";
import styled from "@emotion/styled";
import { PT_TEXT_COLOR } from "src/components/common/theme";

interface Props {
  isSelected: boolean;
}

export const MenuItemWrapper = styled.span<Props>`
  .${Classes.MENU_ITEM} {
    color: ${PT_TEXT_COLOR};
    letter-spacing: -0.1px;
    line-height: 18px;
    padding: 7px 8px;

    .${Classes.FILL} {
      font-weight: ${(props: Props) => (props.isSelected ? 500 : undefined)};
      margin-right: 8px;
    }

    &:focus {
      outline: none;
    }
  }
`;
