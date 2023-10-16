/**
 * Test suite for collection util-related functionality.
 */

// App dependencies
import { expect } from "@playwright/test";
import { getDOIPath } from "src/views/Collection/utils";
import { test } from "tests/common/test";

const { describe } = test;

describe("collection util", () => {
  describe("Get DOI Path", () => {
    test("returns DOI path", () => {
      const doiPath = "123/456";
      const doiUrl = `https://doi.org/${doiPath}`;
      const path = getDOIPath(doiUrl);
      expect(path).toEqual(doiPath);
    });
  });
});
