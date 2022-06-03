import { makeStyles } from "@material-ui/core";
import { scrollbar } from "src/components/common/Filter/common/style";
import { GRAY } from "src/components/common/theme";
import styled from "styled-components";

export const VIEW_LIST_ITEM_HEIGHT = 32;
export const VIEW_LIST_SUBHEADER_HEIGHT = 23;

interface Props {
  panelWidth: number;
}

interface ScrollProps {
  maxHeight: number;
  scrollable: boolean;
}

export const ViewPanel = styled.div<Props>`
  min-width: ${(props) =>
    `${props.panelWidth}px`}; /* required; makes allowances for list item selected state font weight changes by maintaining panel min width */
`;

/* eslint-disable sort-keys -- ignore object key order for style objects */
export const useFilterViewStyles = makeStyles({
  viewDivider: {
    backgroundColor: "rgba(16, 22, 26, 0.15)",
    marginLeft: 8,
    marginRight: 8,
  },
  viewHeading: {
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

export const ViewPanelScroll = styled.div<ScrollProps>`
  max-height: ${({ maxHeight }) => `${maxHeight}px`};
  overflow-y: auto;
  padding-right: ${({ scrollable }) => (scrollable ? "8px" : undefined)};
  ${scrollbar};
`;
