import styled from "@emotion/styled";
import { gray100, grayWhite } from "src/common/theme";
import { CommonThemeProps } from "@czi-sds/components";
import { DRAG_EVENT_TYPE } from "src/components/Collection/components/CollectionDatasetsGrid/components/Row/DatsetRow/components/ReorderModeRow/index";

interface RowProps extends CommonThemeProps {
  dragEventType?: DRAG_EVENT_TYPE;
}

export const Row = styled("tr")<RowProps>`
  cursor: grab;
  margin: ${rowMargin}px;
  padding: ${rowPadding}px;

  &:hover {
    background-color: ${rowBackgroundColor};
  }

  td {
    opacity: ${cellOpacity};
  }
`;

function cellOpacity(props: RowProps): string | undefined {
  if (props.dragEventType === DRAG_EVENT_TYPE.DRAG_START) {
    return "0.8";
  }
  if (props.dragEventType === DRAG_EVENT_TYPE.DRAG) {
    return "0";
  }
}

function rowBackgroundColor(props: RowProps): string | undefined {
  if (props.dragEventType === DRAG_EVENT_TYPE.DRAG) {
    return grayWhite();
  }
  return gray100(props);
}

function rowMargin(props: RowProps): string | undefined {
  if (props.dragEventType === DRAG_EVENT_TYPE.DRAG_START) {
    return "0 -16";
  }
}

function rowPadding(props: RowProps): string | undefined {
  if (props.dragEventType === DRAG_EVENT_TYPE.DRAG_START) {
    return "0 16";
  }
}
