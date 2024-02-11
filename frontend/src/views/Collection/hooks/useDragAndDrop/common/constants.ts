import {
  Dragging,
  DRAGGING_DIRECTION,
  OffsetByIndex,
} from "src/views/Collection/hooks/useDragAndDrop/common/entities";
import { css } from "@emotion/react";

export const DEFAULT_CLIENT_Y = 0;
export const DEFAULT_DIRECTION_TOLERANCE = 8;
export const DEFAULT_DRAGGING: Dragging = {
  draggingDirection: DRAGGING_DIRECTION.DOWN,
  draggingIndex: 0,
  droppingIndex: 0,
  shadowIndex: 0,
  size: 0,
};
export const DEFAULT_DRAGGING_STYLES = css`
  transform: translateY(0);
  transition: none;
`;
export const DEFAULT_OFFSETS: OffsetByIndex = new Map();
export const SELECTOR_TABLE_BODY = "table tbody";
export const SELECTOR_TABLE_ROW = "tr";
