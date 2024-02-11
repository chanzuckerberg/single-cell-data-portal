export interface Dragging {
  draggingDirection: DRAGGING_DIRECTION;
  draggingIndex: number; // The index of the dragged element.
  droppingIndex: number; // The index of the droppable target element.
  shadowIndex: number; // The index of the dragged element that provides a visual cue for the dropping position.
  size: number; // The number of elements.
}

export enum DRAGGING_DIRECTION {
  UP = "UP",
  DOWN = "DOWN",
}

export type OffsetByIndex = Map<number, number>; // Offset keyed by element index.
