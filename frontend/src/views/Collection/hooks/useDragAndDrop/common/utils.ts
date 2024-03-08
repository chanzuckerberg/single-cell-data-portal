import { DragAndDrop } from "src/views/Collection/hooks/useDragAndDrop/common/entities";
import { UseDragAndDrop } from "src/views/Collection/hooks/useDragAndDrop/useDragAndDrop";

/**
 * Returns drag and drop related values and actions.
 * @param dragAndDrop - Drag and drop.
 * @param datasetIndex - Dataset index.
 * @return drag and drop related values and actions.
 */
export function getDragAndDrop(
  dragAndDrop: Omit<UseDragAndDrop, "dragAndDropStyles">,
  datasetIndex: number
): DragAndDrop {
  return {
    ...dragAndDrop,
    datasetIndex,
  };
}
