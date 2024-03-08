import { css, SerializedStyles } from "@emotion/react";
import { useMemo } from "react";
import { Dragging } from "src/views/Collection/hooks/useDragAndDrop/common/entities";

const TRANSITION_DURATION = 0.1; // Transition duration in seconds.

/**
 * Handles the drag and drop styles for the "Drag and Drop" feature for collection view dataset reordering.
 * @returns dragging styles.
 */
export function useDragAndDropStyles(
  dragging: Dragging
): SerializedStyles | undefined {
  const {
    dragTargetIndex,
    dragTargetTranslateY,
    dropTargetIndex,
    dropTargetTranslateY,
  } = dragging;

  return useMemo(
    () =>
      calculateDragAndDropStyles(
        dragTargetIndex,
        dropTargetIndex,
        dragTargetTranslateY,
        dropTargetTranslateY
      ),
    [
      dragTargetIndex,
      dragTargetTranslateY,
      dropTargetIndex,
      dropTargetTranslateY,
    ]
  );
}

/**
 * Returns styles for the dragging and dropping elements during the drag operation.
 * @param dragTargetIndex - Original index of the dragged element in the list.
 * @param dropTargetIndex - Original index of the drop target in the list.
 * @param dragTargetTranslateY - Vertical translation value of the dragged element.
 * @param dropTargetTranslateY - Vertical translation value of the dropped elements.
 * @returns drag and drop styles.
 */
export function calculateDragAndDropStyles(
  dragTargetIndex: number,
  dropTargetIndex: number,
  dragTargetTranslateY?: number,
  dropTargetTranslateY?: number
): SerializedStyles | undefined {
  // Drag operation is not in progress.
  if (!dragTargetTranslateY || !dropTargetTranslateY) return;
  // No change in position - dragged element is in its original position.
  if (dragTargetTranslateY === 0 && dropTargetTranslateY === 0) {
    return getOriginalTransformStyles();
  }
  // Calculate transition styles - dragged element has been moved.
  // Get the directional multiplier for the translation.
  const multiplier = dropTargetIndex > dragTargetIndex ? 1 : -1;
  // Get the drag target selector.
  const dragTargetSelector = dragTargetIndex + 1;
  // Get the drop target selectors.
  const dropTargetSelectors = calculateDropTargetSelectors(
    dragTargetIndex,
    dropTargetIndex
  );
  return getTransformStyles(
    dragTargetSelector,
    dragTargetTranslateY,
    dropTargetSelectors,
    dropTargetTranslateY,
    multiplier
  );
}

/**
 * Returns the drop target selectors, i.e. the nth-of-type selectors of the dropped elements.
 * @param dragTargetIndex - Original index of the dragged element in the list.
 * @param dropTargetIndex - Original index of the drop target in the list.
 * @returns tuple containing the lowest and highest indexes of the dropped elements.
 */
function calculateDropTargetSelectors(
  dragTargetIndex: number,
  dropTargetIndex: number
): [number, number] {
  const dragTargetSelector = dragTargetIndex + 1;
  const dropTargetSelector = dropTargetIndex + 1;
  if (dropTargetSelector < dragTargetSelector) {
    // Drop target is above the original drag target position.
    // Selectors will be from the drop position, to the above the drag target position.
    return [dropTargetSelector, dragTargetSelector - 1];
  } else {
    // Drop target is below the original drag target position.
    // Selectors will be from the below the drag target position, to the drop position.
    return [dragTargetSelector + 1, dropTargetSelector];
  }
}

/**
 * Returns the transform styles for the original position of the dragged and dropped elements.
 * @returns transform style.
 */
function getOriginalTransformStyles() {
  return css`
    transform: translateY(0);
    transition: transform ${TRANSITION_DURATION}s ease-in-out;
  `;
}

/**
 * Returns the transform styles for the dragging and dropping elements.
 * @param dragTargetSelector - Dragged element's nth-of-type selector.
 * @param dragTargetTranslateY - Vertical translation value of the dragged element.
 * @param dropTargetSelectors - Tuple containing the dropped elements' nth-of-type selectors.
 * @param dropTargetTranslateY - Vertical translation value of the dropped elements.
 * @param multiplier - Translation directional multiplier.
 * @returns transform styles.
 */
function getTransformStyles(
  dragTargetSelector: number,
  dragTargetTranslateY: number,
  dropTargetSelectors: [number, number],
  dropTargetTranslateY: number,
  multiplier: number
) {
  const [dropSelectorStart, dropSelectorEnd] = dropTargetSelectors;
  return css`
    transition: transform ${TRANSITION_DURATION}s ease-in-out;

    &:nth-of-type(${dragTargetSelector}) {
      transform: translateY(${multiplier * dragTargetTranslateY}px);
    }

    &:nth-of-type(n + ${dropSelectorStart}):nth-of-type(
        -n + ${dropSelectorEnd}
      ) {
      transform: translateY(${-multiplier * dropTargetTranslateY}px);
    }
  `;
}
