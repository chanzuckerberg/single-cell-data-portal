import { Dragging } from "src/views/Collection/hooks/useDragAndDrop/common/entities";

export const DEFAULT_DIRECTION_TOLERANCE = 8;
export const DEFAULT_DRAGGING: Dragging = {
  dragClientY: 0,
  dragTargetIndex: 0,
  dropTargetIndex: 0,
};
export const SELECTOR_TABLE_BODY = "table tbody";
export const SELECTOR_TABLE_ROW = "tr";
