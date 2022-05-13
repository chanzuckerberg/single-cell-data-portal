import { Classes, Icon } from "@blueprintjs/core";
import styled from "styled-components";

export const Wrapper = styled.div`
  position: relative;
`;

export const StyledDangerIcon = styled(Icon)`
  &.${Classes.ICON} {
    position: absolute;
    right: 8px;
    top: 15px;
  }
`;
