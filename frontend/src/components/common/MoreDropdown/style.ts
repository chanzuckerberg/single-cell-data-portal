import { Classes } from "@blueprintjs/core";
import styled from "@emotion/styled";
import { ActionButton } from "src/components/common/Grid/components/ActionButton/style";
import { GRAY } from "src/components/common/theme";

export const MoreButton = styled(ActionButton)`
  &.${Classes.BUTTON}.${Classes.MINIMAL} svg {
    fill: ${GRAY.A};
  }
`;
