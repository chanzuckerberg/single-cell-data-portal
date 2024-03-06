import styled from "@emotion/styled";
import { gray100, grayWhite, spacesS } from "src/common/theme";
import { CommonThemeProps } from "@czi-sds/components";
import { DRAG_EVENT_TYPE } from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/ReorderModeRow/index";
import { css, SerializedStyles } from "@emotion/react";

interface RowProps extends CommonThemeProps {
  dragEventType?: DRAG_EVENT_TYPE;
}

export const Row = styled("tr")<RowProps>`
  cursor: grab;
  position: relative; /* required; for z-index */
  z-index: 0; /* required; for correct image capture of the row being dragged, when row is positioned under the fixed header */

  &:hover {
    background-color: ${gray100};
  }

  ${draggedRow}
`;

function draggedRow(props: RowProps): SerializedStyles | undefined {
  if (props.dragEventType === DRAG_EVENT_TYPE.DRAG) {
    return css`
      &:hover {
        background-color: ${grayWhite()};
      }

      td {
        opacity: 0;
      }
    `;
  }
  if (props.dragEventType === DRAG_EVENT_TYPE.DRAG_START) {
    return css`
      margin: 0 -${spacesS(props)}px;
      padding: 0 ${spacesS(props)}px;

      td {
        opacity: 0.8;
      }
    `;
  }
}
