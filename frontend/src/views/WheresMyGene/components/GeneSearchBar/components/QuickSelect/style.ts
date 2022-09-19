import styled from "@emotion/styled";
import { Button, ButtonDropdown, Icon, IconButton, MenuItem } from "czifui";


export const ButtonWrapper = styled.div`
  display: flex;
  flex-direction: column;
  gap: 4px;
`;

export const StyledIconButton = styled(Button)`
  min-width: 88px !important;
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