import {
  Dragging,
  DRAGGING_DIRECTION,
} from "src/views/Collection/hooks/useDragAndDrop/common/entities";
import { css } from "@emotion/react";

export const DEFAULT_CLIENT_Y = 0;
export const DEFAULT_DIRECTION = DRAGGING_DIRECTION.DOWN;
export const DEFAULT_DIRECTION_TOLERANCE = 8;
export const DEFAULT_DRAG_AND_DROP_STYLES = css`
  transform: translateY(0);
  transition: none;
`;
export const DEFAULT_DRAGGING: Dragging = {
  dragAndDropIndexes: [],
  dragAndDropStyles: DEFAULT_DRAG_AND_DROP_STYLES,
  draggingIndex: 0,
  droppingIndex: 0,
  offsetByIndex: new Map(),
  shadowIndex: 0,
};
export const SELECTOR_TABLE_BODY = "table tbody";
export const SELECTOR_TABLE_ROW = "tr";
