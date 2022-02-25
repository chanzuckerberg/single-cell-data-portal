import { makeStyles } from "@material-ui/core";
import {
  GRAY,
  LIGHT_GRAY,
  PRIMARY_BLUE,
  PT_TEXT_COLOR,
} from "src/components/common/theme";

/* eslint sort-keys: 0 */ /* TODO(cc) */
export const useSliderStyles = makeStyles({
  root: {
    color: PRIMARY_BLUE,
    display: "block",
    height: 4,
    margin: 0,
    marginLeft: 13, // TODO(cc) temporary
    marginRight: 13, // TODO(cc) temporary
    paddingBottom: 30, // TODO(cc) temporary; use 24
    paddingTop: 34, // TODO(cc) temporary; use 28
    width: 220, // TODO(cc) temporary
    "& .MuiSlider-mark": {
      display: "none",
    },
    "& .MuiSlider-markLabel": {
      color: GRAY.A,
      fontSize: 12,
      letterSpacing: "normal",
      lineHeight: "15px",
      top: 41,
    },
    "& .MuiSlider-rail": {
      borderRadius: 2,
      color: LIGHT_GRAY.C,
      height: 4,
      opacity: 1,
    },
    "& .MuiSlider-thumb": {
      height: 14,
      marginLeft: -7,
      width: 14,
      "&:hover": {
        boxShadow: "none",
      },
    },
    "& .MuiSlider-track": {
      borderRadius: 2,
      height: 4,
    },
    "& .MuiSlider-valueLabel": {
      fontSize: 12,
      fontWeight: 600,
      left: "50%",
      letterSpacing: "normal",
      lineHeight: "15px",
      top: -13,
      transform: "scale(1) translateX(-50%) translateY(-10px) !important",

      "& > span": {
        background: LIGHT_GRAY.E,
        borderRadius: 4,
        height: "unset",
        padding: "2px 4px",
        width: "unset",
      },

      "& *": {
        color: PT_TEXT_COLOR,
        transform:
          "none" /* removes transform rotate from all children of valueLabel */,
      },
    },
  },
});
