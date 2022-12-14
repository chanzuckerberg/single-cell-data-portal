import styled from "@emotion/styled";
import { Button, ButtonIcon, MenuItem } from "czifui";

export const ButtonWrapper = styled.div`
  display: flex;
  flex-direction: column;
`;

export const StyledSelectButton = styled(Button)`
  min-width: 88px !important;
  width: 30px;
  height: 30px;
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
