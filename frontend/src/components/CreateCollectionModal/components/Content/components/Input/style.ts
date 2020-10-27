import { Label } from "@blueprintjs/core";
import { DARK_GRAY } from "src/components/common/theme";
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
