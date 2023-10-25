import {
  DialogActions,
  DialogPaper,
  DialogTitle,
  fontHeaderL,
} from "@czi-sds/components";
import styled from "@emotion/styled";
import { error500, spacesDefault, spacesXl, spacesXxs } from "src/common/theme";
import { Button } from "src/components/common/Button";

export const ActionButton = styled(Button)`
  &:hover {
    background-color: ${error500};
  }
`;

export const StyledDialogTitle = styled(DialogTitle)`
  ${fontHeaderL}
  margin-bottom: ${spacesDefault}px;
`;

export const StyledDialogPaper = styled(DialogPaper)`
  padding: ${spacesXl}px !important;
`;

/**
 * (thuang): We want the word to word space between buttons to be
 * 16px, so reducing margin-left to spacesXxs gives us that
 */
export const StyledDialogAction = styled(DialogActions)`
  &.MuiDialogActions-spacing > :not(:first-of-type) {
    margin-left: ${spacesXxs}px;
  }
`;
