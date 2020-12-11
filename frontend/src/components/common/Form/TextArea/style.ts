import { Icon } from "@blueprintjs/core";
import styled from "styled-components";
import { PT_GRID_SIZE_PX } from "../../theme";

export const Wrapper = styled.div`
  position: relative;
`;

export const StyledIcon = styled(Icon)`
  && {
    position: absolute;
    right: ${PT_GRID_SIZE_PX}px;
    top: ${2 * PT_GRID_SIZE_PX}px;
  }
`;
