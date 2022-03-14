import styled from "@emotion/styled";
import { IconButton, MenuItem } from "czifui";

export const ButtonWrapper = styled.div`
  display: flex;
  flex-direction: column;
  gap: 4px;
  align-items: center;
`;

export const StyledIconButton = styled(IconButton)`
  width: 30px;
  height: 30px;
`;

export const StyledMenuItem = styled(MenuItem)`
  width: 100%;
`;
