import { Classes, Dialog } from "@blueprintjs/core";
import styled from "styled-components";
import { GREY } from "../theme";

export const StyledDialog = styled(Dialog)`
  && {
    .bp3-heading {
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
      background: ${GREY.LIGHT};
      border: 1px solid ${GREY.BORDER_LIGHT};
      border-radius: 4px;
    }

    background-color: white;
    padding-bottom: 0px;
    box-shadow: 0px 2px 4px rgba(0, 0, 0, 0.302065);
    border-radius: 4px;

    /* (thuang): Default Blueprint width */
    min-width: 500px;
    width: unset;
  }
`;
