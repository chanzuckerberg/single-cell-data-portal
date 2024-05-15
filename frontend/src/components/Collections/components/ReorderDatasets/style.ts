import styled from "@emotion/styled";
import {
  spacesS,
  spacesXxs,
  success400,
  success500,
  success600,
} from "src/common/theme";
import { Button } from "src/components/common/Button";

export const SquareButton = styled(Button)`
  background-color: ${success400};
  height: unset;

  &:hover {
    background-color: ${success500};
  }

  &:active {
    background-color: ${success600};
  }
`;

export const MinimalButton = styled(Button)`
  min-width: unset;
  padding: ${spacesXxs}px ${spacesS}px;
`;
