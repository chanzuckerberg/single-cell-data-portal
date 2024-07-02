import styled from "@emotion/styled";
import { Dialog as SDSDialog } from "@czi-sds/components";
import { grey100, shadowL, spacesL } from "src/common/theme";

export const StyledDialog = styled(SDSDialog)`
  .MuiBackdrop-root {
    background-color: rgba(0, 0, 0, 0.3);
  }

  .MuiDialog-paper {
    box-shadow:
      ${shadowL},
      inset 0 0 0 1px ${grey100};
    padding: 32px;
  }

  .MuiDialogTitle-root {
    .MuiTypography-root {
      letter-spacing: -0.46px;
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

export const StyledLink = styled.a`
  width: fit-content;
`;
