import styled from "@emotion/styled";
import { fontCapsXxs, InputText as SDSInputText } from "@czi-sds/components";
import { spacesS, textSecondary } from "src/common/theme";

export const StyledInputText = styled(SDSInputText)`
  &.MuiFormControl-root {
    margin: 0;

    .MuiFormHelperText-root {
      margin: 4px 0 0 0;

      &:first-letter {
        text-transform: capitalize;
      }
    }
  }
`;

export const Label = styled.label`
  ${fontCapsXxs}
  color: ${textSecondary};
  display: block;
  margin-bottom: ${spacesS}px;
`;
