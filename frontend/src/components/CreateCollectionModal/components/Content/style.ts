import { Classes, Divider, H4 } from "@blueprintjs/core";
import { PT_GRID_SIZE_PX, PT_TEXT_COLOR } from "src/components/common/theme";
import styled, { css } from "styled-components";
import Input from "../../../common/Form/Input";

interface Props {
  isFilterEnabled: boolean;
}

export const Form = styled.form<Props>`
  ${(props) => {
    return props.isFilterEnabled
      ? css`
          margin: 0;
          width: 760px;
        `
      : css`
          padding-top: ${PT_GRID_SIZE_PX}px;
          width: 580px;
        `;
  }};
`;

export const CollectionDetail = styled.div`
  display: grid;
  gap: 16px;
  grid-template-columns: 1fr 1fr;
  margin-bottom: 24px;

  /* Collection name and description field */
  [for="description"],
  [for="name"] {
    grid-column: 1 / -1; /* span the grid column specification */
  }
`;

export const CollectionLinks = styled.div`
  display: grid;
  grid-template-columns: 1fr;
  margin: 24px 0 16px;
  row-gap: 16px;
`;

export const Title = styled(H4)`
  color: ${PT_TEXT_COLOR};
  grid-column: 1 / -1; /* span the grid column specification */
  letter-spacing: -0.224px;
  margin: 0;
`;

export const FormDivider = styled(Divider)`
  margin-left: 0;
  margin-right: 0;
`;

/**
 * @deprecated once filter feature flag is removed (#1718).
 */
export const StyledInput = styled(Input)`
  /* Blank for ContactWrapper to target */
`;

/**
 * @deprecated once filter feature flag is removed (#1718).
 */
export const ContactWrapper = styled.div`
  display: flex;

  ${StyledInput}:not(:last-child) {
    margin-right: ${2 * PT_GRID_SIZE_PX}px;
  }
`;

export const CollectionFooter = styled.div`
  margin: 24px 0 0;

  /* Footer actions */
  .${Classes.DIALOG_FOOTER_ACTIONS} {
    gap: 16px;

    /* Footer button */
    .${Classes.BUTTON} {
      margin: 0;
    }
  }
`;
