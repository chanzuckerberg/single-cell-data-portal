import { fontBodyS, fontCapsXxs } from "@czi-sds/components";
import styled from "@emotion/styled";
import { grey500, spacesL, spacesS, spacesXl } from "src/common/theme";
import Loader from "src/components/common/Grid/components/Loader";
import { StyledDialog } from "src/views/Collection/common/style";

export const Dialog = styled(StyledDialog)`
  .MuiDialog-paper {
    gap: ${spacesXl}px;
    min-height: 570px; /* min-height ensures dialog height consistency in any download state. */
  }

  .MuiDialogTitle-root {
    margin: 0;
  }

  .MuiDialogContent-root {
    ${fontBodyS}
    align-content: flex-start;
    display: grid;
    gap: ${spacesXl}px;

    .MuiFormControl-root {
      gap: ${spacesS}px;

      .MuiFormLabel-root {
        ${fontCapsXxs}
        color: ${grey500};
      }

      .MuiFormGroup-row {
        gap: ${spacesL}px;
        margin-bottom: -${spacesL}px; /* negative margin to account for SDS InputRadio FormControlLabel margin see https://github.com/chanzuckerberg/sci-components/blob/main/packages/components/src/core/InputRadio/style.ts#L85 */
      }
    }
  }
`;

export const DialogLoader = styled(Loader)`
  display: flex;

  > div {
    padding: 0;
  }
`;
