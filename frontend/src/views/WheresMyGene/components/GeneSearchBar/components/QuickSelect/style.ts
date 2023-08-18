import styled from "@emotion/styled";
import { autocompleteClasses, menuItemClasses, Popper } from "@mui/material";
import { Button, ButtonIcon, MenuItem } from "@czi-sds/components";
import { OFF_WHITE } from "src/common/theme";

export const ButtonWrapper = styled.div`
  display: flex;
  flex-direction: column;
`;

export const StyledSelectButton = styled(Button)`
  min-width: 88px !important;
  width: 30px;
  height: 30px;
  font-weight: 500;
`;

export const StyledButtonIcon = styled(ButtonIcon)`
  width: 30px;
  height: 30px;
`;

export const StyledMenuItem = styled(MenuItem)`
  width: 100%;
`;

export const StyledButtonText = styled.span`
  margin-right: 5px;
  margin-top: 3px;
`;

export const StyledFixedSizeList = styled.span`
  background-color: pink;
`;

export const StyledPopper = styled(Popper)`
  /* Overwrite default MUI styles */
  && {
    .${autocompleteClasses.focused}.${menuItemClasses.root}[aria-selected="true"] {
      background-color: transparent;
    }

    .${autocompleteClasses.focused}.${menuItemClasses.root}:hover {
      background-color: ${OFF_WHITE};
    }
  }
`;
