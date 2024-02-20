import { UseDragAndDrop } from "src/views/Collection/hooks/useDragAndDrop/useDragAndDrop";

export interface DragAndDrop extends Omit<UseDragAndDrop, "dragging"> {
  datasetIndex: number;
}

export interface Dragging {
  dragClientY: number; // Dragging mouse y-axis position.
  draggingDirection?: DRAGGING_DIRECTION; // The direction of the dragged element.
  dragTargetIndex: number; // Original index of the dragged element in the list.
  dragTargetTranslateY?: number; // Vertical translation value of the dragged element; the sum of the heights of the elements crossed since drag start.
  dropTargetIndex: number; // Original index of the drop target in the list.
  dropTargetTranslateY?: number; // Vertical translation value of the dropped elements; the height of the dragged element.
  elementHeightByIndex?: ElementHeightByIndex; // Element height keyed by element index.
}

export enum DRAGGING_DIRECTION {
  DOWN = "DOWN",
  UP = "UP",
}

export type ElementHeightByIndex = Map<number, number>;
