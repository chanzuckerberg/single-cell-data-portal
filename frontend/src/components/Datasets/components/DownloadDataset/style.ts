import {
  Dialog as SDSDialog,
  fontBodyM,
  fontBodyS,
  fontCapsXxs,
} from "@czi-sds/components";
import styled from "@emotion/styled";
import {
  cornersS,
  fontWeightBold,
  fontWeightSemibold,
  grey100,
  grey500,
  shadowL,
  spacesL,
  spacesM,
  spacesS,
  spacesXl,
  spacesXs,
  spacesXxs,
  spacesXxxs,
  textPrimary,
  textSecondary,
} from "src/common/theme";
import Loader from "src/components/common/Grid/components/Loader";
import { BLUE, OLD_BLUE } from "src/components/common/theme";
import { Classes } from "@blueprintjs/core";

export const Dialog = styled(SDSDialog)`
  .MuiBackdrop-root {
    background-color: rgba(0, 0, 0, 0.3);
  }

  .MuiDialog-paper {
    border: 1px solid ${grey100};
    box-shadow: ${shadowL};
    gap: ${spacesXl}px;
    min-height: 490px;
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

/**
 * @deprecated by Dialog styles once "DOWNLOAD_UX" feature flag is removed.
 */
export const DownloadDialog = styled(SDSDialog)`
  .MuiBackdrop-root {
    background-color: rgba(17, 20, 24, 0.7);
  }

  .MuiDialog-paper {
    box-shadow: 0 ${spacesXxxs}px ${spacesXxs}px rgba(0, 0, 0, 0.3);
    margin: ${spacesXl}px 0;
    min-height: unset;
    padding: 20px 20px 0 20px;
    width: unset;
  }

  .MuiDialogTitle-root {
    align-items: center;
    display: flex;
    margin: 0;
    min-height: 65px;
    padding: ${spacesXxs}px ${spacesM}px;

    .MuiTypography-root {
      ${fontBodyM}
      font-weight: ${fontWeightSemibold};
      line-height: normal;
    }
  }

  .MuiDialogContent-root {
    ${fontBodyS}
    display: flex;
    flex-direction: column;
    gap: 22px;
    height: 300px;
    letter-spacing: normal;
    line-height: inherit; /* line height from body */
    margin: ${spacesM}px;
    overflow: visible;
    width: 625px;

    .MuiFormControl-root {
      display: grid;
      gap: ${spacesS}px;
      min-height: unset;

      .MuiFormLabel-root {
        ${fontCapsXxs}
        color: ${OLD_BLUE};
        font-weight: ${fontWeightBold};
        letter-spacing: normal;
        line-height: 14px;
      }

      .MuiFormGroup-root {
        flex-direction: row;
        gap: ${spacesL}px;
        margin-bottom: ${spacesS}px;

        .MuiFormControlLabel-root {
          display: flex;
          gap: ${spacesS}px;
          margin: 0;
        }

        .MuiFormControlLabel-label {
          font: inherit;
          letter-spacing: normal;

          &.Mui-disabled {
            color: rgba(95, 107, 124, 0.6);
          }
        }

        .MuiRadio-root {
          color: #738091;
          padding: 0;
          position: relative; /* positions pseudo-element styles for checked and disabled */

          &.Mui-checked:after,
          &.Mui-disabled:after {
            border-radius: 50%;
            content: "";
            display: block;
            height: ${spacesL}px;
            position: absolute;
            width: ${spacesL}px;
          }

          &.Mui-checked {
            color: ${OLD_BLUE};

            &:hover {
              color: ${BLUE.B};
            }

            &:after {
              box-shadow: inset 0 0 0 1px rgba(17, 20, 24, 0.2);
            }
          }

          &.Mui-disabled {
            .MuiSvgIcon-root {
              visibility: hidden;
            }

            &:after {
              background-color: rgba(143, 153, 168, 0.15);
              box-shadow: none;
            }
          }

          &.Mui-checked.Mui-disabled {
            color: rgba(45, 114, 210, 0.5);

            .MuiSvgIcon-root {
              visibility: visible;
            }
          }
        }
      }
    }
  }

  .MuiDialogActions-root {
    gap: ${spacesS}px;
    height: 52px;
    margin: 0 ${spacesS}px;

    .MuiButton-root {
      ${fontBodyS}
      border-radius: ${cornersS}px;
      height: 30px;
      letter-spacing: normal;
      line-height: inherit;
      min-width: unset;
      padding: ${spacesXs}px ${spacesM}px;

      &.MuiButton-text {
        color: #1c2127;
        margin: 0;

        &:hover {
          background-color: rgba(143, 153, 168, 0.15);
        }

        &:active {
          background-color: rgba(143, 153, 168, 0.3);
        }
      }

      &.MuiButton-containedPrimary {
        background-color: ${BLUE.C};
        color: white;
        margin: 0;

        &:hover {
          background-color: ${BLUE.B};
        }

        &:active {
          background-color: ${BLUE.A};
        }

        &.Mui-disabled {
          background-color: rgba(14, 125, 236, 0.5);
          color: rgba(255, 255, 255, 0.6);
        }
      }
    }
  }

  .${Classes.SPINNER} {
    justify-content: flex-start;
  }
`;

export const DialogLoader = styled(Loader)`
  display: flex;

  > div {
    padding: 0;
  }
`;
