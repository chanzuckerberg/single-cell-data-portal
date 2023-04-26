import { Classes, Icon } from "@blueprintjs/core";
import { css } from "@emotion/react";
import styled from "@emotion/styled";
import { CommonThemeProps, getColors } from "czifui";
import { PRIMARY_BLUE, PT_TEXT_COLOR } from "src/components/common/theme";

export const Filter = styled.div`
  display: grid;
  gap: 8px;
  margin-bottom: 12px;

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
