export interface Dragging {
  draggingDirection: DRAGGING_DIRECTION;
  draggingIndex: number; // The index of the dragged element.
  droppingIndex: number; // The index of the droppable target element.
  indexes: number[];
  offsetByIndex: OffsetByIndex; // Offsets keyed by element index.
  shadowIndex: number; // The index of the dragged element that provides a visual cue for the dropping position.
}

export enum DRAGGING_DIRECTION {
  UP = "UP",
  DOWN = "DOWN",
}

export type OffsetByIndex = Map<number, number>; // Element offset keyed by element index.
