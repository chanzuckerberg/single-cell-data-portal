import { Dialog } from "@blueprintjs/core";
import styled from "styled-components";

export const StyledDialog = styled(Dialog)`
  && {
    .bp3-heading {
      font-size: 16px;
      line-height: 19px;
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
