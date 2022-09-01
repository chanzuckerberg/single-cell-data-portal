import { Classes, Icon } from "@blueprintjs/core";
import styled from "@emotion/styled";

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
