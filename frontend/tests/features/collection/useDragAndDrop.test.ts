/**
 * Test suite for collection dataset drag and drop hook.
 */

import { test } from "tests/common/test";
import { expect } from "@playwright/test";
import {
  getDraggingDirection,
  updateDragging,
} from "src/views/Collection/hooks/useDragAndDrop/useDragAndDrop";
import { DRAGGING_DIRECTION } from "src/views/Collection/hooks/useDragAndDrop/common/entities";
import { DEFAULT_DRAGGING } from "src/views/Collection/hooks/useDragAndDrop/common/constants";

const { describe } = test;

describe("useDragAndDrop", () => {
  const DRAGGING_INDEX = 3;
  const DRAGGING_INDEXES = [0, 1, 2, 3, 4, 5];
  const DRAGGING_INDEXES_DOWN_ONE = [0, 1, 2, 4, 3, 5];
  const DRAGGING_INDEXES_DOWN_TWO = [0, 1, 2, 4, 5, 3];
  const DRAGGING_INDEXES_UP_ONE = [0, 1, 3, 2, 4, 5];
  const OFFSET_BY_INDEX = new Map([
    [0, 64],
    [1, 64],
    [2, 64],
    [3, 84],
    [4, 84],
    [5, 64],
  ]);
  const DRAGGING = {
    ...DEFAULT_DRAGGING,
    dragAndDropIndexes: DRAGGING_INDEXES,
    draggingIndex: DRAGGING_INDEX,
    droppingIndex: 3,
    offsetByIndex: OFFSET_BY_INDEX,
    shadowIndex: 3,
  };
  const DRAGGING_DOWN_TWO = {
    ...DRAGGING,
    dragAndDropIndexes: DRAGGING_INDEXES_DOWN_TWO,
    droppingIndex: 3,
    shadowIndex: 5,
  };

  describe("Calculate Drag Direction", () => {
    describe("getDraggingDirection", () => {
      test("mouse has not moved", () => {
        const draggingDirection = getDraggingDirection(0);
        expect(draggingDirection).toEqual(undefined);
      });
      test("mouse moved 4 pixels down", () => {
        const draggingDirection = getDraggingDirection(4);
        expect(draggingDirection).toEqual(undefined);
      });
      test("mouse moved 10 pixels down", () => {
        const draggingDirection = getDraggingDirection(10);
        expect(draggingDirection).toEqual(DRAGGING_DIRECTION.DOWN);
      });
      test("mouse moved 10 pixels up", () => {
        const draggingDirection = getDraggingDirection(-10);
        expect(draggingDirection).toEqual(DRAGGING_DIRECTION.UP);
      });
    });
  });

  describe("Update Dragging State", () => {
    describe("updateDragging", () => {
      test("from drag start, dragging down by one position", () => {
        const dragging = updateDragging(DRAGGING, 4, DRAGGING_DIRECTION.DOWN);
        expect(dragging).toEqual(
          expect.objectContaining({
            dragAndDropIndexes: DRAGGING_INDEXES_DOWN_ONE,
            shadowIndex: DRAGGING_INDEXES_DOWN_ONE.indexOf(DRAGGING_INDEX),
          })
        );
      });
      test("from drag start, dragging down by two positions", () => {
        const dragging = updateDragging(DRAGGING, 5, DRAGGING_DIRECTION.DOWN);
        expect(dragging).toEqual(
          expect.objectContaining({
            dragAndDropIndexes: DRAGGING_INDEXES_DOWN_TWO,
            shadowIndex: DRAGGING_INDEXES_DOWN_TWO.indexOf(DRAGGING_INDEX),
          })
        );
      });
      test("from drag start, dragging up by one position", () => {
        const dragging = updateDragging(DRAGGING, 2, DRAGGING_DIRECTION.UP);
        expect(dragging).toEqual(
          expect.objectContaining({
            dragAndDropIndexes: DRAGGING_INDEXES_UP_ONE,
            shadowIndex: DRAGGING_INDEXES_UP_ONE.indexOf(DRAGGING_INDEX),
          })
        );
      });
      test("from two positions down, drag up by one position", () => {
        const dragging = updateDragging(
          DRAGGING_DOWN_TWO,
          5,
          DRAGGING_DIRECTION.UP
        );
        expect(dragging).toEqual(
          expect.objectContaining({
            dragAndDropIndexes: DRAGGING_INDEXES_DOWN_ONE,
            shadowIndex: DRAGGING_INDEXES_DOWN_ONE.indexOf(DRAGGING_INDEX),
          })
        );
      });
    });
  });
});
