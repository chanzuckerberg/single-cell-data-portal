/**
 * Test suite for collection dataset reorder hook.
 */

import { test } from "tests/common/test";
import { expect } from "@playwright/test";
import { buildOrderedIds } from "src/views/Collection/hooks/useReorder/useReorder";

const { describe } = test;

describe("useReorder", () => {
  const DATASETS = [
    "f09db59f-c2cd-4da0-ac40-272ef03ec915",
    "9fdb4989-7153-4f24-a71f-3a9edb30c746",
    "5c02fad1-6ca6-436c-8761-d6c1f268e9ea",
    "8322fe3e-d12c-4abf-8670-4f9fc91bd362",
    "29a461c0-2269-413e-8822-2af87e31dcc8",
    "b10db3ce-2d31-4eec-ba25-b090df32edd4",
  ];
  const DATASETS_ORDERED = [
    "f09db59f-c2cd-4da0-ac40-272ef03ec915",
    "9fdb4989-7153-4f24-a71f-3a9edb30c746",
    "8322fe3e-d12c-4abf-8670-4f9fc91bd362",
    "29a461c0-2269-413e-8822-2af87e31dcc8",
    "5c02fad1-6ca6-436c-8761-d6c1f268e9ea",
    "b10db3ce-2d31-4eec-ba25-b090df32edd4",
  ];
  describe("Update Dataset Order", () => {
    describe("buildOrderedIds", () => {
      test("undefined dataset IDs", () => {
        const orderedIds = buildOrderedIds(undefined, 2, 4);
        expect(orderedIds).toBeUndefined();
      });
      test("dataset ID and target ID are equal", () => {
        const orderedIds = buildOrderedIds(DATASETS, 2, 2);
        expect(orderedIds).toEqual(DATASETS);
      });
      test("dataset ID and target ID are not equal", () => {
        const orderedIds = buildOrderedIds(DATASETS, 2, 4);
        expect(orderedIds).toEqual(DATASETS_ORDERED);
      });
    });
  });
});
