import { useCallback, useState } from "react";
import {
  Dragging,
  DRAGGING_DIRECTION,
  ElementHeightByIndex,
} from "src/views/Collection/hooks/useDragAndDrop/common/entities";
import { OnReorderFn } from "src/views/Collection/hooks/useReorder/useReorder";
import {
  DEFAULT_DIRECTION_TOLERANCE,
  DEFAULT_DRAGGING,
} from "src/views/Collection/hooks/useDragAndDrop/common/constants";
import { HEADER_ID } from "src/components/Header/constants";
import { SerializedStyles } from "@emotion/react";
import { useDragAndDropStyles } from "src/views/Collection/hooks/useDragAndDropStyles";

export interface DragAndDropAction {
  onDraggingOver: (dropTargetIndex: number, clientY: number) => void;
  onDropping: (onReorder: OnReorderFn) => void;
  onEndDragging: () => void;
  onStartDragging: (
    dragTargetIndex: number,
    clientY: number,
    elementHeightByIndex: ElementHeightByIndex
  ) => void;
}

export interface UseDragAndDrop {
  dragAndDropAction: DragAndDropAction;
  dragAndDropStyles?: SerializedStyles;
}

/**
 * "Drag and Drop" feature for collection view dataset reordering.
 * Handles the drag and drop UI/UX for reordering datasets in the collection view.
 * @returns dragging actions and styles.
 */
export function useDragAndDrop(): UseDragAndDrop {
  const [dragging, setDragging] = useState<Dragging>(DEFAULT_DRAGGING);
  const { dragTargetIndex, dropTargetIndex } = dragging;
  const dragAndDropStyles = useDragAndDropStyles(dragging);

  // Callback fired when dragged element is being dragged over a valid drop target.
  const onDraggingOver = useCallback(
    (dropTargetIndex: number, clientY: number) => {
      setDragging((dragging) =>
        calculateNextDraggingState(dragging, dropTargetIndex, clientY)
      );
    },
    []
  );

  // Callback fired when dragged element is dropped on a valid drop target.
  const onDropping = useCallback(
    (onReorder: OnReorderFn) => {
      removeHeaderStyle();
      setDragging(DEFAULT_DRAGGING);
      onReorder(dragTargetIndex, dropTargetIndex);
    },
    [dragTargetIndex, dropTargetIndex]
  );

  // Callback fired when the drag operation ends.
  const onEndDragging = useCallback(() => {
    removeHeaderStyle();
    setDragging(DEFAULT_DRAGGING);
  }, []);

  // Callback fired when dragging starts.
  const onStartDragging = useCallback(
    (
      dragTargetIndex: number,
      clientY: number,
      elementHeightByIndex: ElementHeightByIndex
    ) => {
      setHeaderStyle();
      setDragging((dragging) =>
        calculateStartDraggingState(
          dragging,
          dragTargetIndex,
          clientY,
          elementHeightByIndex
        )
      );
    },
    []
  );

  return {
    dragAndDropAction: {
      onDraggingOver,
      onDropping,
      onEndDragging,
      onStartDragging,
    },
    dragAndDropStyles,
  };
}

/**
 * Updates the dragging mouse y-axis position, based on the change in the y-axis position between consecutive onDraggingOver
 * events. If the change in the position is not significant enough, the last known position is retained.
 * @param draggingState - Current state values relating to drag.
 * @param clientY - Drag event y-axis position.
 * @returns dragging mouse y-axis position.
 */
function calculateDragClientY(draggingState: Dragging, clientY: number) {
  const dy = clientY - draggingState.dragClientY;
  if (shouldCalculateDirection(dy)) {
    return clientY;
  }
  return draggingState.dragClientY;
}

/**
 * Calculates the dragging direction, based on the change in the dragging mouse y-axis positions between consecutive
 * onDraggingOver events. If the change in the position is not significant enough to determine a direction, the last known
 * direction is returned.
 * @param draggingState - Current state values relating to drag.
 * @param clientY - Drag event y-axis position.
 * @returns dragging direction.
 */
export function calculateEffectiveDraggingDirection(
  draggingState: Dragging,
  clientY: number
): DRAGGING_DIRECTION | undefined {
  const dy = clientY - draggingState.dragClientY;
  // Check if the movement is significant enough to calculate a new direction.
  if (shouldCalculateDirection(dy)) {
    // Determine the direction based on the sign of the change in y position.
    return Math.sign(dy) > 0 ? DRAGGING_DIRECTION.DOWN : DRAGGING_DIRECTION.UP;
  }
  // Return the last known direction if the movement is not significant.
  return draggingState.draggingDirection;
}

/**
 * Calculate the next dragging state with changes to the drop target index, and the dragging direction.
 * @param draggingState - Current state values relating to drag.
 * @param dropTargetIndex - Original index of the drop target in the list.
 * @param clientY - Drag event y-axis position.
 * @returns next dragging state.
 */
export function calculateNextDraggingState(
  draggingState: Dragging,
  dropTargetIndex: number,
  clientY: number
): Dragging {
  if (draggingState.dragTargetIndex === dropTargetIndex) {
    // Dragging element is on the drop target; transitioning is complete.
    return draggingState;
  }
  // Calculate the dragging direction.
  const draggingDirection = calculateEffectiveDraggingDirection(
    draggingState,
    clientY
  );
  // Calculate the dragging mouse y-axis position.
  const dragClientY = calculateDragClientY(draggingState, clientY);
  // Adjust the drop target index, based on direction and position of drop target index, relative to the drag target index.
  const adjustedDropTargetIndex = getDropTargetIndex(
    draggingState,
    dropTargetIndex,
    draggingDirection
  );
  // Calculate the y-axis translation values for the drop target elements.
  const dropTargetTranslateY = getDropTargetTranslateYValue(
    draggingState,
    adjustedDropTargetIndex
  );
  // Calculate the y-axis translation values for the drag target elements.
  const dragTargetTranslateY = getDragTargetTranslateYValue(
    draggingState,
    adjustedDropTargetIndex
  );
  return {
    ...draggingState,
    dragClientY,
    draggingDirection,
    dragTargetTranslateY,
    dropTargetIndex: adjustedDropTargetIndex,
    dropTargetTranslateY,
  };
}

/**
 * Calculate the start dragging state.
 * @param draggingState - Current state values relating to drag.
 * @param dragTargetIndex - Original index of the dragged element in the list.
 * @param clientY - Drag event y-axis position.
 * @param elementHeightByIndex - Element height by index.
 * @returns start dragging state.
 */
function calculateStartDraggingState(
  draggingState: Dragging,
  dragTargetIndex: number,
  clientY: number,
  elementHeightByIndex: ElementHeightByIndex
): Dragging {
  return {
    ...draggingState,
    dragClientY: clientY,
    dragTargetIndex,
    dragTargetTranslateY: 0, // Neutral position.
    dropTargetIndex: dragTargetIndex, // Drop target is currently on the drag target.
    dropTargetTranslateY: 0, // Neutral position.
    elementHeightByIndex,
  };
}

/**
 * Returns the drag target y-axis translation value; facilitates the visual positioning of the drag target during dragging action.
 * The drag target is transformed vertically by the sum of the heights of the elements crossed since drag start.
 * @param draggingState - Current state values relating to drag.
 * @param dropTargetIndex - Original index of the drop target in the list.
 * @returns drag target y-axis translation values.
 */
function getDragTargetTranslateYValue(
  draggingState: Dragging,
  dropTargetIndex: number
): number | undefined {
  const { dragTargetIndex, elementHeightByIndex } = draggingState;
  if (!elementHeightByIndex) return; // Drag operation has not started.
  if (dragTargetIndex === dropTargetIndex) return 0; // Dragging element is on the drop target; no transformation needed.
  const startIndex = Math.min(dragTargetIndex, dropTargetIndex);
  const endIndex = Math.max(dragTargetIndex, dropTargetIndex);
  let translateY = 0;
  for (let i = startIndex; i <= endIndex; i++) {
    if (i === dragTargetIndex) continue; // Exclude the dragging element's height.
    translateY += elementHeightByIndex.get(i) ?? 0;
  }
  return translateY;
}

/**
 * Returns the drop target index, adjusted, if required, when direction and position of drop target index relative to
 * the drag target index are considered.
 * Ensures visual consistency with the direction of the drag and intended drop position.
 * @param draggingState - Current state values relating to drag.
 * @param dropTargetIndex - Original index of the drop target in the list.
 * @param draggingDirection - Dragging direction.
 * @returns drop target index.
 */
function getDropTargetIndex(
  draggingState: Dragging,
  dropTargetIndex: number,
  draggingDirection?: DRAGGING_DIRECTION
): number {
  if (
    isDraggingUpwardsWithDropTargetBelow(
      draggingState,
      dropTargetIndex,
      draggingDirection
    )
  ) {
    // The direction is upwards and the drop target is currently positioned below the original drag target position.
    // The drop target position should be adjusted to above the given drop target index.
    return dropTargetIndex - 1;
  }
  if (
    isDraggingDownwardsWithDropTargetAbove(
      draggingState,
      dropTargetIndex,
      draggingDirection
    )
  ) {
    // The direction is downwards and the drop target is currently positioned above the original drag target position.
    // The drop target position should be adjusted to below the given drop target index.
    return dropTargetIndex + 1;
  }
  // Drop target position is consistent with the direction of the drag and drop position and no adjustment is required.
  // i.e. direction is downwards and the drop target is currently positioned below the original drag target position, or
  // direction is upwards and the drop target is currently positioned above the original drag target position.
  return dropTargetIndex;
}

/**
 * Returns the drop target y-axis translation value; facilitates the visual positioning of the drop targets during dragging action.
 * The drop targets are transformed vertically by the height of the dragged element.
 * @param draggingState - Current state values relating to drag.
 * @param dropTargetIndex - Original index of the drop target in the list.
 * @returns drop target y-axis translation values.
 */
function getDropTargetTranslateYValue(
  draggingState: Dragging,
  dropTargetIndex: number
): number | undefined {
  const { dragTargetIndex, elementHeightByIndex } = draggingState;
  if (!elementHeightByIndex) return; // Drag operation has not started.
  if (dragTargetIndex === dropTargetIndex) return 0; // Dragging element is on the drop target; no transformation needed.
  return elementHeightByIndex.get(dragTargetIndex) || 0;
}

/**
 * Returns true if the dragging direction is downwards and the drop target is currently positioned above the original
 * drag target.
 * @param draggingState - Current state values relating to drag.
 * @param dropTargetIndex - Original index of the drop target in the list.
 * @param draggingDirection - Dragging direction.
 * @returns true if the dragging direction is downwards and the drop target is currently positioned above the original drag target.
 */
function isDraggingDownwardsWithDropTargetAbove(
  draggingState: Dragging,
  dropTargetIndex: number,
  draggingDirection?: DRAGGING_DIRECTION
) {
  return (
    draggingDirection === DRAGGING_DIRECTION.DOWN &&
    dropTargetIndex < draggingState.dragTargetIndex
  );
}

/**
 * Returns true if the dragging direction is upwards and the drop target is currently positioned below the original
 * drag target.
 * @param draggingState - Current state values relating to drag.
 * @param dropTargetIndex - Original index of the drop target in the list.
 * @param draggingDirection - Dragging direction.
 * @returns true if the dragging direction is upwards and the drop target is currently positioned below the original drag target.
 */
function isDraggingUpwardsWithDropTargetBelow(
  draggingState: Dragging,
  dropTargetIndex: number,
  draggingDirection?: DRAGGING_DIRECTION
) {
  return (
    draggingDirection === DRAGGING_DIRECTION.UP &&
    dropTargetIndex > draggingState.dragTargetIndex
  );
}

/**
 * Removes the header style to enable pointer events.
 */
function removeHeaderStyle() {
  const headerEl = document.getElementById(HEADER_ID);
  if (headerEl) {
    headerEl.style.removeProperty("pointer-events");
  }
}

/**
 * Sets the header style to disable pointer events.
 * Facilitates the upward scroll action, while dragging.
 */
function setHeaderStyle() {
  const headerEl = document.getElementById(HEADER_ID);
  if (headerEl) {
    headerEl.style.setProperty("pointer-events", "none");
  }
}

/**
 * Determines whether the movement along the y-axis is sufficient to warrant a direction calculation.
 * This decision is based on comparing the magnitude of the change in y-axis position (delta y) against a defined
 * tolerance and ensures that minor movements do not trigger direction changes and / or calculations.
 * @param dy - Delta dragging y-axis position.
 * @param tolerance - Tolerance for detecting the direction change.
 * @returns true if the change in the y-axis position is sufficient to detect a direction change.
 */
function shouldCalculateDirection(
  dy: number,
  tolerance = DEFAULT_DIRECTION_TOLERANCE
) {
  return Math.abs(dy) > tolerance;
}
