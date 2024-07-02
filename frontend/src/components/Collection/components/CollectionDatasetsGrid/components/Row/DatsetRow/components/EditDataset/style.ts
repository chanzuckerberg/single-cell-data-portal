import styled from "@emotion/styled";
import { StyledDialog as CommonStyledDialog } from "src/views/Collection/common/style";

export const StyledDialog = styled(CommonStyledDialog)`
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
