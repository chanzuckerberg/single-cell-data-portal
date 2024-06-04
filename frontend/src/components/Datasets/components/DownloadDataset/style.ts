import {
  Dialog as SDSDialog,
  fontBodyS,
  fontCapsXxs,
} from "@czi-sds/components";
import styled from "@emotion/styled";
import {
  grey100,
  grey500,
  shadowL,
  spacesL,
  spacesS,
  spacesXl,
} from "src/common/theme";
import Loader from "src/components/common/Grid/components/Loader";

export const Dialog = styled(SDSDialog)`
  .MuiBackdrop-root {
    background-color: rgba(0, 0, 0, 0.3);
  }

  .MuiDialog-paper {
    border: 1px solid ${grey100};
    box-shadow: ${shadowL};
    gap: ${spacesXl}px;
    min-height: 570px; /* min-height ensures dialog height consistency in any download state. */
    padding: 32px;
  }

  .MuiDialogTitle-root {
    margin: 0;

    .MuiTypography-root {
      letter-spacing: -0.46px;
    }
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

  .MuiDialogActions-root {
    gap: ${spacesL}px;
    margin: 0;
    padding-top: ${spacesL}px;

    .MuiButton-root {
      &.MuiButton-containedPrimary {
        margin: 0;
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
