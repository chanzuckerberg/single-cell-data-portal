import { Classes, Divider, H4 } from "@blueprintjs/core";
import styled from "@emotion/styled";
import { PT_TEXT_COLOR } from "src/components/common/theme";
import { DropdownForm } from "src/components/common/Form/Dropdown/style";

export const Form = styled.form`
  margin: 0;
  width: 760px;
`;

export const CollectionDetail = styled.div`
  display: grid;
  gap: 16px;
  grid-template-columns: 1fr 1fr;
  margin-bottom: 24px;

  /* Collection name and description field, and consortia dropdown */
  [for="description"],
  [for="name"],
  ${DropdownForm} {
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
