import { Icon, InputGroup, Label } from "@blueprintjs/core";
import { DARK_GRAY, PT_GRID_SIZE_PX } from "src/components/common/theme";
import styled from "styled-components";

interface StyledLProps {
  percentage?: number;
}

/**
 * @deprecated - supersede by StyledFormLabel once filter feature flag is removed (#1718).
 */
export const StyledLabel = styled(Label)`
  ${({ percentage = 100 }: StyledLProps) => {
    return `width: ${percentage}%`;
  }}
`;

/**
 * @deprecated - supersede by FormLabelText once filter feature flag is removed (#1718).
 */
export const LabelText = styled.span`
  color: ${DARK_GRAY.C};
`;

/**
 * @deprecated - supersede by StyledFormLabel once filter feature flag is removed (#1718).
 */
export const StyledIcon = styled(Icon)`
  margin-right: ${PT_GRID_SIZE_PX}px;
`;

/**
 * @deprecated - supersede by StyledFormLabel once filter feature flag is removed (#1718).
 */
export const StyledInputGroup = styled(InputGroup)`
  && input {
    border-radius: unset;
  }
`;

/**
 * @deprecated - supersede by SelectFormLabel once filter feature flat is removed (#1718).
 */
export const StyledDiv = styled.div`
  display: flex;
  flex-direction: column;
`;
