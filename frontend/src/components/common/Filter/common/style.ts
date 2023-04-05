import { Classes, Icon } from "@blueprintjs/core";
import { css } from "@emotion/react";
import styled from "@emotion/styled";
import { CommonThemeProps, getColors } from "czifui";
import { PRIMARY_BLUE, PT_TEXT_COLOR } from "src/components/common/theme";

const gray300 = (props: CommonThemeProps) => getColors(props)?.gray[300];

export const Filter = styled.div`
  display: grid;
  font-feature-settings: normal; /* required; overrides layout.css specification */
  gap: 4px;
  margin-bottom: 8px;

  &:last-child {
    margin-bottom: 0;
  }

  /* Filter open with "active" button. */

  & .${Classes.POPOVER_OPEN} {
    .${Classes.BUTTON} {
      color: ${PT_TEXT_COLOR};
    }
  }
`;

export const FilterDivider = styled.hr`
  background: none;
  box-shadow: inset 0px -0.5px 0px ${gray300};
  height: 0.5px;
  margin: 16px 0;
`;

export const SelectionIcon = styled(Icon)`
  && {
    align-items: center;
    color: ${PRIMARY_BLUE};
    display: flex;
    height: 18px;
    justify-content: center;
    margin: 0 8px 0 0;

    > svg {
      height: 14px;
      width: 14px;
    }
  }
`;

export const scrollbar = (props: CommonThemeProps) => {
  const colors = getColors(props);
  return css`
    &::-webkit-scrollbar {
      width: 4px;
    }

    &::-webkit-scrollbar-thumb {
      background-clip: content-box;
      background-color: ${colors?.gray["200"]};
      border-radius: 4px;

      &:hover {
        background-color: ${colors?.gray["300"]};
      }
    }
  `;
};

export const InfoButtonWrapper = styled.span`
  padding-left: 4px;
  cursor: pointer;
`;

export const NotificationWrapper = styled.div`
  position: absolute;
  z-index: 999;
  overflow: hidden;
`;
