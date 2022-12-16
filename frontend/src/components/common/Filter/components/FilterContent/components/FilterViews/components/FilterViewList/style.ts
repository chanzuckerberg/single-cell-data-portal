import { css } from "@emotion/react";
import styled from "@emotion/styled";
import { ListItemText } from "@mui/material";
import { List, ListItem, ListItemProps } from "czifui";
import { GRAY } from "src/components/common/theme";

interface Props {
  selected: boolean;
}

export const ViewSublist = styled(List)`
  margin-left: 22px;
`;

// (thuang): MUIv5 deprecated `button` prop, but not removed yet
// We should add ListItemButton to SDS and use that instead
// https://mui.com/material-ui/api/list-item/#props
interface TempListItemProps extends ListItemProps {
  button?: boolean;
}

const StyledListItem = styled(ListItem)<TempListItemProps>`
  /* TODO(cc) remove && after updating SDS version that has this commit https://github.com/chanzuckerberg/sci-components/pull/201 */
  && {
    letter-spacing: -0.1px;
    line-height: 18px;
    margin: 0 /* overrides margin from layout.css */;
    padding: 7px 8px;

    &:before {
      display: none; /* remove list item bullet. */
    }
  }
`;

export const NoMatches = styled(StyledListItem)`
  /* TODO(cc) remove && after updating SDS version that has this commit https://github.com/chanzuckerberg/sci-components/pull/201 */
  && {
    color: ${GRAY.A};
  }
`;

export const ViewListItem = styled(StyledListItem)`
  /* TODO(cc) remove && after updating SDS version that has this commit https://github.com/chanzuckerberg/sci-components/pull/201 */
  && {
    &:hover {
      background-color: rgba(
        167,
        182,
        194,
        0.3
      ); /* mimics bp menu item pseudo-class hover background color */
    }
    &:active {
      background-color: rgba(
        115,
        134,
        148,
        0.3
      ); /* mimics bp menu item pseudo-class active background color */
    }
  }
`;

export const ViewListItemText = styled(ListItemText)<Props>`
  display: flex;
  margin: 0;

  /* List item text - "primary" */
  span:first-of-type {
    flex: 1;
    ${({ selected }) => {
      return (
        selected &&
        css`
          font-weight: 500;
        `
      );
    }}
    margin-right: 8px;
  }

  /* List item count - "secondary" */
  span:last-of-type {
    color: ${GRAY.A};
  }
`;
