/**
 * Test suite for collection dataset drag and drop hook.
 */

import { test } from "tests/common/test";
import { expect } from "@playwright/test";
import {
  calculateEffectiveDraggingDirection,
  calculateNextDraggingState,
} from "src/views/Collection/hooks/useDragAndDrop/useDragAndDrop";
import {
  Dragging,
  DRAGGING_DIRECTION,
} from "src/views/Collection/hooks/useDragAndDrop/common/entities";
import { DEFAULT_DRAGGING } from "src/views/Collection/hooks/useDragAndDrop/common/constants";

const { describe } = test;

describe("useDragAndDrop", () => {
  const DRAG_CLIENT_Y_START = 232; // Assume mid-point of the dragged element, and top of first element is 0.
  const DRAG_TARGET_INDEX = 3;
  const DROP_TARGET_INDEX = 3;
  const ELEMENT_HEIGHT_BY_INDEX = new Map([
    [0, 64],
    [1, 64],
    [2, 64],
    [3, 84],
    [4, 84],
    [5, 64],
  ]);
  const DRAGGING_STATE: Dragging = {
    ...DEFAULT_DRAGGING,
    dragClientY: DRAG_CLIENT_Y_START,
    dragTargetIndex: DRAG_TARGET_INDEX,
    dropTargetIndex: DROP_TARGET_INDEX,
    elementHeightByIndex: ELEMENT_HEIGHT_BY_INDEX,
  };

  describe("Calculate Effective Drag Direction", () => {
    describe("calculateEffectiveDraggingDirection", () => {
      test("mouse has not moved", () => {
        const draggingDirection = calculateEffectiveDraggingDirection(
          DRAGGING_STATE,
          DRAG_CLIENT_Y_START
        );
        expect(draggingDirection).toBeUndefined();
      });
      test("mouse moved 4 pixels down", () => {
        const draggingDirection = calculateEffectiveDraggingDirection(
          DRAGGING_STATE,
          DRAG_CLIENT_Y_START + 4
        );
        expect(draggingDirection).toBeUndefined();
      });
      test("mouse moved 10 pixels down", () => {
        const draggingDirection = calculateEffectiveDraggingDirection(
          DRAGGING_STATE,
          DRAG_CLIENT_Y_START + 10
        );
        expect(draggingDirection).toEqual(DRAGGING_DIRECTION.DOWN);
      });
      test("mouse moved 10 pixels up", () => {
        const draggingDirection = calculateEffectiveDraggingDirection(
          DRAGGING_STATE,
          DRAG_CLIENT_Y_START - 10
        );
        expect(draggingDirection).toEqual(DRAGGING_DIRECTION.UP);
      });
    });
  });

  describe("Calculate Next Dragging State", () => {
    describe("calculateNextDraggingState", () => {
      test("from drag start, drag down by one position", () => {
        const dropTargetIndex = DROP_TARGET_INDEX + 1;
        const dragging = calculateNextDraggingState(
          DRAGGING_STATE,
          dropTargetIndex,
          276
        );
        expect(dragging).toEqual(
          expect.objectContaining({
            draggingDirection: DRAGGING_DIRECTION.DOWN,
            dragTargetTranslateY: 84,
            dropTargetIndex: dropTargetIndex,
            dropTargetTranslateY: 84,
          })
        );
      });
      test("from drag start, drag down by two positions", () => {
        const dropTargetIndex = DROP_TARGET_INDEX + 2;
        const dragging = calculateNextDraggingState(
          DRAGGING_STATE,
          dropTargetIndex,
          360
        );
        expect(dragging).toEqual(
          expect.objectContaining({
            draggingDirection: DRAGGING_DIRECTION.DOWN,
            dragTargetTranslateY: 148,
            dropTargetIndex: dropTargetIndex,
            dropTargetTranslateY: 84,
          })
        );
      });
      test("from drag start, drag up by one position", () => {
        const dropTargetIndex = DROP_TARGET_INDEX - 1;
        const dragging = calculateNextDraggingState(
          DRAGGING_STATE,
          dropTargetIndex,
          192
        );
        expect(dragging).toEqual(
          expect.objectContaining({
            draggingDirection: DRAGGING_DIRECTION.UP,
            dragTargetTranslateY: 64,
            dropTargetIndex: dropTargetIndex,
            dropTargetTranslateY: 84,
          })
        );
      });
      test("from two positions down, drag up by one position", () => {
        const CURRENT_DRAGGING_STATE = {
          ...DRAGGING_STATE,
          draggingDirection: DRAGGING_DIRECTION.DOWN,
          dragClientY: 392,
          dropTargetIndex: 5,
        };
        const dragging = calculateNextDraggingState(
          CURRENT_DRAGGING_STATE,
          5,
          360
        );
        expect(dragging).toEqual(
          expect.objectContaining({
            draggingDirection: DRAGGING_DIRECTION.UP,
            dragTargetTranslateY: 84,
            dropTargetIndex: 4,
            dropTargetTranslateY: 84,
          })
        );
      });
    });
  });
});
