import { Classes, Icon } from "@blueprintjs/core";
import styled from "styled-components";
import { PT_GRID_SIZE_PX } from "../../theme";

export const Wrapper = styled.div`
  position: relative;
`;

/**
 * @deprecated - supersede by StyledDangerIcon once filter feature flag is removed (#1718).
 */
export const StyledIcon = styled(Icon)`
  && {
    position: absolute;
    right: ${PT_GRID_SIZE_PX}px;
    top: ${2 * PT_GRID_SIZE_PX}px;
  }
`;

export const StyledDangerIcon = styled(Icon)`
  &.${Classes.ICON} {
    position: absolute;
    right: 8px;
    top: 8px;
  }
`;
