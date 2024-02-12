import { DragEvent, useCallback, useRef, useState } from "react";
import {
  Dragging,
  DRAGGING_DIRECTION,
  OffsetByIndex,
} from "src/views/Collection/hooks/useDragAndDrop/common/entities";
import { css, SerializedStyles } from "@emotion/react";
import { OnReorderFn } from "src/views/Collection/hooks/useReorderMode";
import {
  DEFAULT_CLIENT_Y,
  DEFAULT_DIRECTION,
  DEFAULT_DIRECTION_TOLERANCE,
  DEFAULT_DRAGGING,
} from "src/views/Collection/hooks/useDragAndDrop/common/constants";

export interface DragAndDropAction {
  onDragging: (dragEvent: DragEvent<HTMLElement>) => void;
  onDraggingOver: (droppingIndex: number) => void;
  onDropping: (onReorder: OnReorderFn) => void;
  onEndDragging: () => void;
  onStartDragging: (
    draggingIndex: number,
    offsetByIndex: OffsetByIndex
  ) => void;
}

interface UseDragAndDrop {
  dragAndDropAction: DragAndDropAction;
  dragAndDropStyles?: SerializedStyles;
}

/**
 * "Drag and Drop" feature for collection view dataset reordering.
 * Handles the drag and drop UI/UX for reordering datasets in the collection view.
 * @returns drag and drop actions and dragging styles.
 */
export function useDragAndDrop(): UseDragAndDrop {
  const clientYRef = useRef<number>(DEFAULT_CLIENT_Y);
  const draggingDirectionRef = useRef<DRAGGING_DIRECTION>(DEFAULT_DIRECTION);
  const [dragging, setDragging] = useState<Dragging>(DEFAULT_DRAGGING);
  const { dragAndDropStyles, draggingIndex, shadowIndex } = dragging;

  // Callback fired every few milliseconds when the dragged element is being dragged.
  // Dragging direction is updated based on the change in the y-axis position.
  // Direction is used to position the visual cue of the dragged element.
  const onDragging = useCallback((dragEvent: DragEvent<HTMLElement>) => {
    const clientY = dragEvent.clientY;
    // Update dragging direction.
    const direction = getDraggingDirection(clientY - clientYRef.current);
    if (!direction) return;
    // Update clientY position, only when direction has been detected.
    clientYRef.current = clientY;
    draggingDirectionRef.current = direction;
  }, []);

  // Callback fired when dragged element is being dragged over a valid drop target.
  const onDraggingOver = useCallback((droppingIndex: number) => {
    setDragging((dragging) =>
      updateDragging(dragging, droppingIndex, draggingDirectionRef.current)
    );
  }, []);

  // Callback fired when dragged element is dropped on a valid drop target.
  const onDropping = useCallback(
    (onReorder: OnReorderFn) => {
      setDragging(DEFAULT_DRAGGING);
      onReorder(draggingIndex, shadowIndex);
    },
    [draggingIndex, shadowIndex]
  );

  // Callback fired when the drag operation ends.
  const onEndDragging = useCallback(() => {
    clientYRef.current = DEFAULT_CLIENT_Y;
    draggingDirectionRef.current = DEFAULT_DIRECTION;
    setDragging(DEFAULT_DRAGGING);
  }, []);

  // Callback fired when dragging starts.
  const onStartDragging = useCallback(
    (draggingIndex: number, offsetByIndex: OffsetByIndex) => {
      setDragging((dragging) =>
        updateDraggingStart(dragging, draggingIndex, offsetByIndex)
      );
    },
    []
  );

  return {
    dragAndDropAction: {
      onDragging,
      onDraggingOver,
      onDropping,
      onEndDragging,
      onStartDragging,
    },
    dragAndDropStyles,
  };
}

/**
 * Returns the dragging and dropping y-offset values.
 * @param draggingIndex - Dragging index.
 * @param indexes - Drag indexes affected by the drag operation.
 * @param offsetByIndex - Offset, keyed by index.
 * @returns dragging and dropping y-offset values.
 */
function getDragAndDropOffset(
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
  const droppingOffset =
    draggingOffset === 0 ? 0 : offsetByIndex.get(draggingIndex) ?? 0;
  return [draggingOffset, droppingOffset];
}

/**
 * Updates the dragging and dropping styles for the dragged element and eligible drop elements.
 * @param draggingIndex - Dragging index.
 * @param shadowIndex - Shadow index.
 * @param offsetByIndex - Offset by index.
 * @returns drag and drop styles.
 */
function getDragAndDropStyles(
  draggingIndex: number,
  shadowIndex: number,
  offsetByIndex: OffsetByIndex
): SerializedStyles {
  const draggingOffsetDirection = shadowIndex > draggingIndex ? 1 : -1;
  const indexes = getDraggingIndexes(draggingIndex, shadowIndex);
  const offset = getDragAndDropOffset(draggingIndex, indexes, offsetByIndex);
  // Transform dragging and dropping elements on the y-axis.
  const [draggingOffset, droppingOffset] = offset;
  // No change in position - shadowIndex is within the bounds of the original dragging index.
  if (draggingOffset === 0 && droppingOffset === 0) {
    return css`
      transform: translateY(0);
      transition: transform 0.1s ease-in-out;
    `;
  }
  const [nthStart, nthEnd] = getSelectorPositions(draggingIndex, indexes);
  return css`
    transition: transform 0.1s ease-in-out;

    &:nth-of-type(${draggingIndex + 1}) {
      transform: translateY(${draggingOffset * draggingOffsetDirection}px);
    }

    &:nth-of-type(n + ${nthStart}):nth-of-type(-n + ${nthEnd}) {
      transform: translateY(${droppingOffset * -draggingOffsetDirection}px);
    }
  `;
}

/**
 * Returns the dragging direction.
 * @param dy - Delta clientY.
 * @param tolerance - Tolerance for detecting the direction change.
 * @returns direction.
 */
function getDraggingDirection(
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
 * Returns the dropping index, adjusted, if required, by direction.
 * A downward direction should result in the dropping index being the maximum value between the dropping index, and
 * an upward direction should result in the dropping index being the minimum value.
 * @param draggingState - Dragging state.
 * @param droppingIndex - Dropping index.
 * @param draggingDirection - Dragging direction.
 * @returns droppingIndex.
 */
function getDroppingIndex(
  draggingState: Dragging,
  droppingIndex: number,
  draggingDirection: DRAGGING_DIRECTION
): number {
  if (draggingDirection === DRAGGING_DIRECTION.UP) {
    return Math.min(droppingIndex, draggingState.droppingIndex);
  }
  return Math.max(droppingIndex, draggingState.droppingIndex);
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
 * Reorders the drag and drop indexes with a revised dropping index.
 * @param draggingIndex - Dragging index.
 * @param droppingIndex - Dropping index.
 * @param dragAndDropIndexes - Drag and drop indexes.
 * @returns drag and drop indexes.
 */
function reorderDragAndDropIndexes(
  draggingIndex: number,
  droppingIndex: number,
  dragAndDropIndexes: number[]
): number[] {
  const updatedDragAndDropIndexes = [...dragAndDropIndexes];
  const shadowIndex = dragAndDropIndexes.indexOf(draggingIndex);
  // Remove the shadow index (i.e. dragging index) from its current position.
  updatedDragAndDropIndexes.splice(shadowIndex, 1);
  // Insert the dragging index to the dropping index.
  updatedDragAndDropIndexes.splice(droppingIndex, 0, draggingIndex);
  return updatedDragAndDropIndexes;
}

/**
 * Updates the dragging state with changes to the dropping index, and the dragging direction.
 * @param draggingState - Dragging state.
 * @param droppingIndex - Dropping index.
 * @param draggingDirection - Dragging direction.
 * @returns updated dragging state.
 */
function updateDragging(
  draggingState: Dragging,
  droppingIndex: number,
  draggingDirection: DRAGGING_DIRECTION
): Dragging {
  if (draggingState.draggingIndex === droppingIndex) {
    // Dragging element is on the drop target; transitioning is complete.
    return { ...draggingState, droppingIndex, isTransitioning: false };
  }
  if (draggingState.isTransitioning) {
    // During the transition of the dragging element, the synchronization between the shadow index and the dropping index
    // should be maintained to ensure consistency in the visual representation of the drag operation.
    // Aggressive dragging actions, during transitions, can cause discrepancies between the two indexes.
    // By setting the isTransitioning flag to false, we signal the completion of this specific transition, prompting a
    // recalibration of the shadow index.
    if (draggingState.shadowIndex !== draggingState.droppingIndex) {
      return { ...draggingState, isTransitioning: false };
    }
    // During the transition of the dragging element, specifically when the shadow index overlaps with the dropping index,
    // the dropping index value must be maintained. As the shadow index moves over the dropping index, the dropping index
    // will dynamically adjust its position to reflect the transitioning state. However, this adjustment leads to
    // unwanted consequences; the dropping index's index position no longer represents the original given value, and
    // dragging styles during this transition result in visual glitches as the transition oscillates
    // between the original dropping index and the dynamically updated dropping index. To avoid these issues, while
    // the transition is ongoing, we effectively "lock" the dropping index to its original value.
    return {
      ...draggingState,
      droppingIndex: getDroppingIndex(
        draggingState,
        droppingIndex,
        draggingDirection
      ),
    };
  }
  // Dragging over a new dropping index.
  const droppingIndexIndex =
    draggingState.dragAndDropIndexes.indexOf(droppingIndex);
  // Reorder the drag and drop indexes with the updated dropping index.
  const dragAndDropIndexes = reorderDragAndDropIndexes(
    draggingState.draggingIndex,
    droppingIndexIndex,
    draggingState.dragAndDropIndexes
  );
  // Shadow index to provide a visual cue for the dropping position.
  const shadowIndex = droppingIndexIndex;
  // Update the drag and drop styles.
  const dragAndDropStyles = getDragAndDropStyles(
    draggingState.draggingIndex,
    shadowIndex,
    draggingState.offsetByIndex
  );
  return {
    ...draggingState,
    dragAndDropIndexes,
    dragAndDropStyles,
    droppingIndex: droppingIndexIndex,
    shadowIndex,
    isTransitioning: true,
  };
}

/**
 * Returns start dragging state.
 * @param draggingState - Dragging state.
 * @param draggingIndex - Dragging index.
 * @param offsetByIndex - Offset by index.
 * @returns dragging state.
 */
function updateDraggingStart(
  draggingState: Dragging,
  draggingIndex: number,
  offsetByIndex: OffsetByIndex
): Dragging {
  return {
    ...draggingState,
    dragAndDropIndexes: [...offsetByIndex.keys()],
    draggingIndex,
    droppingIndex: draggingIndex,
    offsetByIndex,
    shadowIndex: draggingIndex,
  };
}
