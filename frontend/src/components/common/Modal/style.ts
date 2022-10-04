import { Classes, Dialog } from "@blueprintjs/core";
import styled from "@emotion/styled";

export const StyledDialog = styled(Dialog)`
  && {
    .bp4-heading {
      font-size: 16px;
      height: 57px;
      display: flex;
      align-items: center;
      border-radius: 4px;
    }

    .${Classes.DIALOG_FOOTER} {
      display: flex;
      align-items: center;
      justify-content: flex-end;

      margin: 0;
      height: 52px;
    }

    background-color: white;
    padding-bottom: 0px;
    box-shadow: 0px 2px 4px rgba(0, 0, 0, 0.302065);
    border-radius: 4px;

    /* (thuang): Default Blueprint width */
    min-width: 500px;
    width: unset;
  }
  /* Dataset download loading state */
  &.modal-loading {
    min-height: 461px;
    min-width: 697px;
  }
`;
