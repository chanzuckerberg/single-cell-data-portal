import {
  GRAY,
  LIGHT_GRAY,
  PRIMARY_BLUE,
  PT_TEXT_COLOR,
} from "src/components/common/theme";
import Slider from "@mui/material/Slider";
import styled from "@emotion/styled";

export const StyledSlider = styled(Slider)`
  color: ${PRIMARY_BLUE};
  display: block;
  height: 4px;
  margin: 0;
  margin-left: 24px;
  margin-right: 24px;
  padding-bottom: 30px;
  padding-top: 34px;
  width: 220px;

  & .MuiSlider-mark {
    display: none;
  }

  & .MuiSlider-markLabel {
    color: ${GRAY.A};
    font-size: 12px;
    letter-spacing: normal;
    line-height: 15px;
    top: 47px;
  }

  & .MuiSlider-rail {
    border-radius: 2px;
    color: ${LIGHT_GRAY.C};
    height: 4px;
    opacity: 1;
  }

  & .MuiSlider-thumb {
    height: 14px;
    margin-left: -7px;
    width: 14px;

    &:focus {
      outline: none;
    }

    &:hover {
      box-shadow: none;
    }
  }

  & .MuiSlider-track {
    border-radius: 2px;
    height: 4px;
  }

  & .MuiSlider-valueLabel {
    font-size: 12px;
    font-weight: 500;
    left: 50%;
    letter-spacing: normal;
    line-height: 15px;
    top: -13px;
    transform: scale(1) translateX(-50%) translateY(-10px) !important;

    & > span {
      background: ${LIGHT_GRAY.E};
      border-radius: 4px;
      height: unset;
      padding: 2px 4px;
      width: unset;
    }

    & * {
      color: ${PT_TEXT_COLOR};
      transform: none; /* removes transform rotate from all children of valueLabel */
    }
  }
`;
