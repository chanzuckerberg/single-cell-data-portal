import styled from "@emotion/styled";
import { ListItemText } from "@material-ui/core";
import { List, ListItem } from "czifui";
import { GRAY } from "src/components/common/theme";

export const ViewSublist = styled(List)`
  margin-left: 22px;
`;

const StyledListItem = styled(ListItem)`
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

export const ViewListItemText = styled(ListItemText)`
  display: flex;
  margin: 0;
`;
