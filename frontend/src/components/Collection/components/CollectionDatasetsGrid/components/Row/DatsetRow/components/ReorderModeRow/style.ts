import styled from "@emotion/styled";
import { gray100, grayWhite } from "src/common/theme";
import { CommonThemeProps } from "@czi-sds/components";
import { css } from "@emotion/react";
import { DRAG_EVENT_TYPE } from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/ReorderModeRow/index";

interface RowProps extends CommonThemeProps {
  dragEventType?: DRAG_EVENT_TYPE;
}

export const Row = styled("tr")<RowProps>`
  cursor: grab;

  &:active {
    cursor: grabbing;
  }

  &:hover {
    background-color: ${gray100};
  }

  ${(props) => {
    return (
      props.dragEventType === DRAG_EVENT_TYPE.DRAG_START &&
      css`
        & {
          margin: 0 -16px;
          padding: 0 16px;

          td {
            opacity: 0.8;
          }
        }
      `
    );
  }};

  ${(props) => {
    return (
      props.dragEventType === DRAG_EVENT_TYPE.DRAG &&
      css`
        & {
          &:hover {
            background-color: ${grayWhite()};
          }

          td {
            opacity: 0;
          }
        }
      `
    );
  }};
`;
