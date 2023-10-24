import {
  Button,
  DialogTitle,
  fontBodyS,
  fontHeaderL,
} from "@czi-sds/components";
import styled from "@emotion/styled";
import { error500, spacesDefault, textSecondary } from "src/common/theme";

const FONT_WEIGHT_500 = 500;

export const ActionButton = styled(Button)`
  ${fontBodyS}
  font-weight: ${FONT_WEIGHT_500};

  &:hover {
    background-color: ${error500};
  }
`;

export const CancelButton = styled(Button)`
  ${fontBodyS}
  color: ${textSecondary};
  font-weight: ${FONT_WEIGHT_500};
`;

export const StyledDialogTitle = styled(DialogTitle)`
  ${fontHeaderL}
  margin-bottom: ${spacesDefault}px;
`;
