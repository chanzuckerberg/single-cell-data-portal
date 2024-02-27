import styled from "@emotion/styled";
import { Menu as SDSMenu, MenuItem as SDSMenuItem } from "@czi-sds/components";
import {
  error400,
  grey400,
  shadowM,
  spacesS,
  spacesXxs,
} from "src/common/theme";

export const Menu = styled(SDSMenu)`
  &.MuiPopover-root {
    z-index: 10; /* positions menu beneath BP portals to allow proper layering with delete/publish collection alert interactions */
  }

  .MuiMenu-paper {
    box-shadow: ${shadowM};
    min-width: 199px; /* matches mocks */
  }
`;

export const MenuItem = styled(SDSMenuItem)`
  &.MuiMenuItem-root {
    .primary-text {
      gap: ${spacesS}px;

      > span {
        display: flex;
        margin: ${spacesXxs}px 0;
      }
    }
  }
`;

export const DeleteMenuItem = styled(MenuItem)`
  &.MuiMenuItem-root {
    .primary-text {
      color: ${error400};
    }
  }
`;

export const ReorderMenuItem = styled(MenuItem)`
  &.MuiMenuItem-root {
    svg {
      color: ${grey400};
      margin: ${spacesXxs}px 0;
    }

    &.Mui-disabled {
      svg {
        color: inherit;
      }
    }
  }
`;
