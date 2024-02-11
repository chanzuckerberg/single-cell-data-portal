import { DragEvent, useCallback, useEffect, useRef, useState } from "react";
import {
  Dragging,
  DRAGGING_DIRECTION,
  OffsetByIndex,
} from "src/views/Collection/hooks/useDragAndDrop/common/entities";
import { css, SerializedStyles } from "@emotion/react";
import { OnReorderFn } from "src/views/Collection/hooks/useReorderMode";
import {
  DEFAULT_CLIENT_Y,
  DEFAULT_DIRECTION_TOLERANCE,
  DEFAULT_DRAGGING,
  DEFAULT_DRAGGING_STYLES,
  DEFAULT_OFFSETS,
  SELECTOR_TABLE_BODY,
  SELECTOR_TABLE_ROW,
} from "src/views/Collection/hooks/useDragAndDrop/common/constants";

export interface DragAndDropAction {
  onDragging: (dragEvent: DragEvent<HTMLElement>) => void;
  onDraggingOver: (droppingIndex: number) => void;
  onDropping: (onReorder: OnReorderFn) => void;
  onEndDragging: () => void;
  onStartDragging: (
    dragEvent: DragEvent<HTMLElement>,
    draggingIndex: number
  ) => void;
}

interface UseDragAndDrop {
  dragAndDropAction: DragAndDropAction;
  draggingStyles?: SerializedStyles;
}

/**
 * "Drag and Drop" feature for collection view dataset reordering.
 * Handles the drag and drop UI/UX for reordering datasets in the collection view.
 * @returns drag and drop actions and dragging styles.
 */
export function useDragAndDrop(): UseDragAndDrop {
  const clientYRef = useRef<number>(DEFAULT_CLIENT_Y);
  const offsetByIndexRef = useRef<OffsetByIndex>(DEFAULT_OFFSETS);
  const [dragging, setDragging] = useState<Dragging>(DEFAULT_DRAGGING);
  const [draggingStyles, setDraggingStyles] = useState<SerializedStyles>();
  const { draggingIndex, droppingIndex, shadowIndex } = dragging;

  // Callback fired every few milliseconds when the dragged element is being dragged.
  const onDragging = useCallback((dragEvent: DragEvent<HTMLElement>) => {
    // Update dragging direction.
    const clientY = dragEvent.clientY;
    const direction = calculateDirection(clientY - clientYRef.current);
    if (!direction) return;
    clientYRef.current = clientY; // Update clientY reference only when direction has been detected.
    setDragging((dragging) => updateDraggingDirection(dragging, direction));
  }, []);

  // Callback fired when dragged element is being dragged over a valid drop target.
  const onDraggingOver = useCallback((droppingIndex: number) => {
    setDragging((dragging) =>
      updateDraggingDroppingIndex(dragging, droppingIndex)
    );
  }, []);

  // Callback fired when dragged element is dropped on a valid drop target.
  const onDropping = useCallback(
    (onReorder: OnReorderFn) => {
      setDraggingStyles(DEFAULT_DRAGGING_STYLES);
      onReorder(draggingIndex, droppingIndex);
    },
    [draggingIndex, droppingIndex]
  );

  // Callback fired when the drag operation ends.
  const onEndDragging = useCallback(() => {
    clientYRef.current = DEFAULT_CLIENT_Y;
    offsetByIndexRef.current = DEFAULT_OFFSETS;
    setDragging(DEFAULT_DRAGGING);
  }, []);

  // Callback fired when dragging starts.
  const onStartDragging = useCallback(
    (dragEvent: DragEvent<HTMLElement>, draggingIndex: number) => {
      offsetByIndexRef.current = getOffsets(dragEvent);
      setDragging((dragging) => ({
        ...dragging,
        draggingIndex,
        droppingIndex: draggingIndex,
        shadowIndex: draggingIndex,
        size: offsetByIndexRef.current.size,
      }));
    },
    []
  );

  // Update dragging state with changes in dropping index.
  useEffect(() => {
    setDragging(updateDragging);
  }, [droppingIndex]);

  // Update the dragging styles with changes to the shadow index.
  useEffect(() => {
    setDraggingStyles(
      updateDraggingStyles(draggingIndex, shadowIndex, offsetByIndexRef.current)
    );
  }, [draggingIndex, shadowIndex]);

  return {
    dragAndDropAction: {
      onDragging,
      onDraggingOver,
      onDropping,
      onEndDragging,
      onStartDragging,
    },
    draggingStyles: draggingStyles,
  };
}

/**
 * Returns the drag direction.
 * @param dy - Delta position.
 * @param tolerance - Tolerance for detecting the direction change.
 * @returns direction.
 */
function calculateDirection(
  dy: number,
  tolerance = DEFAULT_DIRECTION_TOLERANCE
): DRAGGING_DIRECTION | undefined {
  if (dy === 0) return;
  const absY = Math.abs(dy);
  if (absY < tolerance) return; // Tolerance threshold for detecting direction change.
  const sign = Math.sign(dy);
  return sign > 0 ? DRAGGING_DIRECTION.DOWN : DRAGGING_DIRECTION.UP;
}

/**
 * Returns the dragging and dropping y-offset values.
 * @param draggingIndex - Dragging index.
 * @param indexes - Drag indexes affected by the drag operation.
 * @param offsetByIndex - Offset, keyed by index.
 * @returns dragging and dropping y-offset values.
 */
function calculateOffset(
  draggingIndex: number,
  indexes: number[],
  offsetByIndex: OffsetByIndex
): [number, number] {
  // Dragging element transforms on the y-axis by the displaced dropping elements.
  let draggingOffset = 0;
  for (const i of indexes) {
    if (i === draggingIndex) continue;
    draggingOffset += offsetByIndex.get(i) ?? 0;
  }
  // Dropping targets transform on the y-axis by the displaced dragging element.
  const dropOffset =
    draggingOffset === 0 ? 0 : offsetByIndex.get(draggingIndex) ?? 0;
  return [draggingOffset, dropOffset];
}

/**
 * Returns the shadow index (the shadow provides visual positioning of the dragged element during dragging).
 * @param draggingState - Dragging state.
 * @returns shadow index.
 */
function calculateShadowIndex(draggingState: Dragging): number {
  const { draggingDirection, droppingIndex, shadowIndex, size } = draggingState;
  let nextShadowIndex = shadowIndex;
  if (draggingDirection === DRAGGING_DIRECTION.DOWN) {
    if (droppingIndex > shadowIndex) {
      nextShadowIndex = droppingIndex;
    }
  } else {
    if (droppingIndex < shadowIndex) {
      nextShadowIndex = droppingIndex;
    }
  }
  return normalizeShadowIndex(nextShadowIndex, size);
}

/**
 * Returns the indexes included in the drag operation.
 * @param draggingIndex - Dragging index.
 * @param shadowIndex - Shadow index.
 * @returns indexes included in the drag operation.
 */
function getDraggingIndexes(
  draggingIndex: number,
  shadowIndex: number
): number[] {
  const minIndex = Math.min(draggingIndex, shadowIndex);
  const maxIndex = Math.max(draggingIndex, shadowIndex);
  if (minIndex === maxIndex) return [minIndex];
  return Array.from(
    { length: maxIndex - minIndex + 1 },
    (_, index) => minIndex + index
  );
}

/**
 * Returns a map of element height by index.
 * Referenced by dragging styles state to calculate the y-axis translation of the dragging element and dropping elements
 * during drag operation.
 * @param dragEvent - Drag event.
 * @returns a map of offset by index.
 */
function getOffsets(dragEvent: DragEvent<HTMLElement>): OffsetByIndex {
  const tbodyEl = dragEvent.currentTarget.closest(SELECTOR_TABLE_BODY);
  const rowsEl = tbodyEl?.querySelectorAll(SELECTOR_TABLE_ROW);
  if (!rowsEl) return DEFAULT_OFFSETS;
  const offsetByIndex = new Map();
  rowsEl.forEach((rowEl, index) => {
    offsetByIndex.set(index, rowEl.getBoundingClientRect().height);
  });
  return offsetByIndex;
}

/**
 * Returns the start and end position of the nth selectors for the dropping elements.
 * @param draggingIndex - Dragging index.
 * @param indexes - Drag indexes affected by the drag operation.
 * @returns selector positions.
 */
function getSelectorPositions(
  draggingIndex: number,
  indexes: number[]
): [number, number] {
  const nArray = [...indexes];
  nArray.splice(indexes.indexOf(draggingIndex), 1);
  const firstIndex = nArray.shift() || 0;
  const lastIndex = nArray.pop() || firstIndex;
  const start = firstIndex + 1;
  const end = lastIndex + 1;
  return [start, end];
}

/**
 * Normalizes shadow index to remain within the allowable range of the indexes.
 * @param shadowIndex - Shadow index.
 * @param size - Total number of elements.
 * @returns normalized shadow index.
 */
function normalizeShadowIndex(shadowIndex: number, size: number): number {
  if (shadowIndex < 0) return 0;
  const maxIndex = size - 1;
  if (shadowIndex > maxIndex) return maxIndex;
  return shadowIndex;
}

/**
 * Updates the dragging state with the recalculated position of the shadow index.
 * @param draggingState - Dragging state.
 * @returns updated dragging state.
 */
function updateDragging(draggingState: Dragging): Dragging {
  if (draggingState.size === 0) return draggingState; // Dragging is ended and dragging state is reset.
  const shadowIndex = calculateShadowIndex(draggingState);
  return { ...draggingState, shadowIndex };
}

/**
 * Updates the dragging styles.
 * @param draggingIndex - Dragging index.
 * @param shadowIndex - Shadow index.
 * @param offsetByIndex - Offset by index.
 * @returns dragging styles.
 */
function updateDraggingStyles(
  draggingIndex: number,
  shadowIndex: number,
  offsetByIndex: OffsetByIndex
): SerializedStyles {
  const draggingOffsetDirection = shadowIndex > draggingIndex ? 1 : -1;
  const indexes = getDraggingIndexes(draggingIndex, shadowIndex);
  const offset = calculateOffset(draggingIndex, indexes, offsetByIndex);
  // Transform dragging and dropping elements on the y-axis.
  const [draggingOffset, droppingOffset] = offset;
  // No change in position - shadowIndex is within the bounds of the original dragging index.
  if (draggingOffset === 0 && droppingOffset === 0) {
    return css`
      transform: translateY(0);
      transition: transform 0.2s ease-in-out;
    `;
  }
  const [nthStart, nthEnd] = getSelectorPositions(draggingIndex, indexes);
  return css`
    transition: transform 0.2s ease-in-out;

    &:nth-of-type(${draggingIndex + 1}) {
      transform: translateY(${draggingOffset * draggingOffsetDirection}px);
    }

    &:nth-of-type(n + ${nthStart}):nth-of-type(-n + ${nthEnd}) {
      transform: translateY(${droppingOffset * -draggingOffsetDirection}px);
    }
  `;
}

/**
 * Updates the dragging state with changes to the dragging direction.
 * @param draggingState - Dragging state.
 * @param direction - Direction.
 * @returns updated dragging state.
 */
function updateDraggingDirection(
  draggingState: Dragging,
  direction: DRAGGING_DIRECTION
): Dragging {
  if (draggingState.draggingDirection === direction) return draggingState;
  return { ...draggingState, draggingDirection: direction };
}

/**
 * Updates the dragging state with changes to the dropping index.
 * @param draggingState - Dragging state.
 * @param droppingIndex - Dropping index.
 * @returns updated dragging state.
 */
function updateDraggingDroppingIndex(
  draggingState: Dragging,
  droppingIndex: number
): Dragging {
  if (draggingState.draggingIndex === droppingIndex) {
    return draggingState; // Dragging element is the drop target (usually onDragStart or after transition).
  }
  return {
    ...draggingState,
    droppingIndex: droppingIndex,
  };
}
