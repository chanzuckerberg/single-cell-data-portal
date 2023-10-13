/**
 * Test suite for select filter-related utils.
 */

import { expect } from "@playwright/test";
import { buildAnalyticsPayload } from "src/common/hooks/useCategoryFilter/common/selectUtils";
import { test } from "tests/common/test";

const { describe } = test;

describe("filter", () => {
  const KEY_PAYLOAD = "payload";
  const KEY_TISSUE = "tissue";
  const VALUE_TISSUE = "brain";
  describe("Analytics", () => {
    test("builds analytics payload with no payload key", () => {
      const payload = buildAnalyticsPayload(VALUE_TISSUE);
      expect(payload[KEY_PAYLOAD]).toBeTruthy();
      expect(payload[KEY_PAYLOAD]).toEqual(VALUE_TISSUE);
    });
    test("builds analytics payload with payload key", () => {
      const payload = buildAnalyticsPayload(VALUE_TISSUE, KEY_TISSUE);
      expect(payload[KEY_TISSUE]).toBeTruthy();
      expect(payload[KEY_TISSUE]).toEqual(VALUE_TISSUE);
    });
  });
});
