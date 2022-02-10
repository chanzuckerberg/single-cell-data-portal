import { Classes } from "@blueprintjs/core";
import { ActionButton } from "src/components/common/Grid/components/ActionButton/style";
import { GRAY } from "src/components/common/theme";
import styled from "styled-components";

export const MoreButton = styled(ActionButton)`
  &.${Classes.BUTTON}.${Classes.MINIMAL} svg {
    fill: ${GRAY.A};
  }
`;
