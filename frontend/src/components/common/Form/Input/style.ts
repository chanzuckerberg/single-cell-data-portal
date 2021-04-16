import { Icon, InputGroup, Label } from "@blueprintjs/core";
import { DARK_GRAY, PT_GRID_SIZE_PX } from "src/components/common/theme";
import styled from "styled-components";

interface StyledLProps {
  percentage?: number;
}

export const StyledLabel = styled(Label)`
  ${({ percentage = 100 }: StyledLProps) => {
    return `width: ${percentage}%`;
  }}
`;

export const LabelText = styled.span`
  color: ${DARK_GRAY.C};
`;

export const IconWrapper = styled.div`
  display: flex;
  align-items: center;
  height: 100%;
`;

export const StyledIcon = styled(Icon)`
  margin-right: ${PT_GRID_SIZE_PX}px;
`;

export const StyledInputGroup = styled(InputGroup)`
  && input {
    border-radius: unset;
  }
`;

export const StyledDiv = styled.div`
  display: flex;
  flex-direction: column;
`;
