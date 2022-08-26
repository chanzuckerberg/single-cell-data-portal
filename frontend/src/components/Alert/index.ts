import { Alert, Classes } from "@blueprintjs/core";
import styled from "@emotion/styled";
import { DARK_GRAY, PT_GRID_SIZE_PX } from "src/components/common/theme";

export default styled(Alert)`
  width: ${PT_GRID_SIZE_PX * 55}px;
  line-height: 18px;
  .${Classes.BUTTON} {
    :not(.${Classes.INTENT_DANGER}):not(.${Classes.INTENT_PRIMARY}) {
      background: none;
      box-shadow: none !important;
      color: ${DARK_GRAY.A};

      &:hover {
        background: rgba(167, 182, 194, 0.3);
      }
    }
  }
`;
