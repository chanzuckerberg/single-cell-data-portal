import styled from "@emotion/styled";
import { StyledDialog } from "src/views/Collection/common/style";

export const Dialog = styled(StyledDialog)`
  .MuiDialog-container {
    .MuiPaper-root {
      display: block;
      min-height: unset;

      .MuiDialogContent-root {
        margin: 24px 0;
      }
    }
  }
`;
