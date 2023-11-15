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
  spacesM,
  spacesS,
  spacesXl,
  spacesXs,
  textPrimary,
  textSecondary,
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

      .MuiFormGroup-root {
        flex-direction: row;
        gap: ${spacesL}px;

        .MuiFormControlLabel-root {
          display: flex;
          gap: ${spacesS}px;
          margin: 0;
        }

        .MuiRadio-root {
          padding: 0;
        }
      }
    }
  }

  .MuiDialogActions-root {
    gap: ${spacesL}px;
    margin: 0;
    padding-top: ${spacesL}px;

    .MuiButton-root {
      ${fontBodyS}
      font-weight: 500;
      height: unset;
      min-width: unset;

      &.MuiButton-text {
        color: ${textSecondary};

        &:hover,
        &:active {
          color: ${textPrimary};
        }
      }

      &.MuiButton-containedPrimary {
        margin: 0;
        padding: ${spacesXs}px ${spacesM}px;
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
