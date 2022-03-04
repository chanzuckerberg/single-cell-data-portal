import { makeStyles } from "@material-ui/core";
import { PT_TEXT_COLOR } from "src/components/common/theme";

/* eslint-disable sort-keys -- ignore object key order for style objects */
export const useFilterPanelListStyles = makeStyles({
  listItem: {
    color: PT_TEXT_COLOR,
    letterSpacing: "-0.1px",
    lineHeight: "18px",
    margin: 0 /* overrides margin from layout.css */,
    padding: "7px 8px",
    "&:hover": {
      backgroundColor: "rgba(167, 182, 194, 0.3)",
      /* mimics bp menu item pseudo-class hover background color */
    },
    "&:active": {
      backgroundColor: "rgba(115, 134, 148, 0.3)",
      /* mimics bp menu item pseudo-class active background color */
    },
  },
  listItemText: {
    display: "flex",
    margin: 0,
  },
  sublist: {
    marginLeft: 22,
  },
});
/* eslint-enable sort-keys -- ignore object key order for style objects */
