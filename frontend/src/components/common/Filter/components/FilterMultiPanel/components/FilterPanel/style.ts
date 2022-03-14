import { makeStyles } from "@material-ui/core";
import { scrollbar } from "src/components/common/Filter/common/style";
import { GRAY } from "src/components/common/theme";
import styled from "styled-components";

const LIST_HEIGHT = 32;
const LIST_SUBHEADER_HEIGHT = 23;
export const MAX_DISPLAYABLE_LIST_ITEMS = 15;
const PANEL_MAX_HEIGHT =
  (MAX_DISPLAYABLE_LIST_ITEMS + 0.5) * LIST_HEIGHT + LIST_SUBHEADER_HEIGHT;

interface Props {
  panelWidth: number;
  scrollable: boolean;
}

export const Panel = styled.div<Props>`
  max-height: ${PANEL_MAX_HEIGHT}px;
  min-width: ${(props) =>
    `${props.panelWidth}px`}; /* required; makes allowances for list item selected state font weight changes by maintaining panel min width */
  overflow-y: auto;
  padding-right: ${(props) => (props.scrollable ? "8px" : undefined)};
  ${scrollbar};
`;

/* eslint-disable sort-keys -- ignore object key order for style objects */
export const useFilterPanelStyles = makeStyles({
  panelDivider: {
    backgroundColor: "rgba(16, 22, 26, 0.15)",
    marginLeft: 8,
    marginRight: 8,
  },
  panelHeading: {
    color: GRAY.A,
    cursor: "default",
    fontSize: 11,
    fontWeight: 500,
    letterSpacing: "0.03em",
    lineHeight: "15px",
    marginBottom: 8,
    paddingLeft: 8,
    textTransform: "uppercase",
  },
});
/* eslint-enable sort-keys -- ignore object key order for style objects */
