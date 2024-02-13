import { SerializedStyles } from "@emotion/react";

export interface Dragging {
  dragAndDropIndexes: number[];
  dragAndDropStyles: SerializedStyles;
  draggingIndex: number; // The index of the dragged element.
  droppingIndex: number; // The index of the droppable target element.
  offsetByIndex: OffsetByIndex; // Offsets keyed by element index.
  shadowIndex: number; // The index of the dragged element that provides a visual cue for the dropping position.
}

export enum DRAGGING_DIRECTION {
  DOWN = "DOWN",
  UP = "UP",
}

export type OffsetByIndex = Map<number, number>; // Element offset keyed by element index.
