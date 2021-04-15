import { Button, Classes } from "@blueprintjs/core";
import Input from "src/components/common/Form/Input";
import { PT_GRID_SIZE_PX } from "src/components/common/theme";
import styled from "styled-components";

export const Wrapper = styled.div`
  display: flex;
`;

export const IconWrapper = styled.div`
  position: relative;
  width: ${3.5 * PT_GRID_SIZE_PX}px;
`;

export const StyledButton = styled(Button)`
  && {
    position: absolute;
    top: 20px;
    right: -${PT_GRID_SIZE_PX}px;
  }
`;

export const StyledLinkTypeButton = styled(Button)`
  width: 100%;
  height: 34px;
  margin-top: 4px;
  justify-content: space-between;
`;

export const StyledURLInput = styled(Input)``;
export const StyledNameInput = styled(Input)``;

export const LinkWrapper = styled.div`
  display: flex;

  & > * :not(:last-child) {
    margin-right: ${2 * PT_GRID_SIZE_PX}px;
  }

  & > .${Classes.POPOVER_WRAPPER} {
    width: 25%;
  }

  & .${Classes.POPOVER_TARGET} {
    width: 100%;
  }
`;
