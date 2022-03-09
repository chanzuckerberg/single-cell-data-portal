/**
 * Test suite for collection util-related functionality.
 */

// App dependencies
import { getDOIPath } from "src/views/Collection/utils";

describe("collection util", () => {
  describe("Get DOI Path", () => {
    it("returns DOI path", () => {
      const doiPath = "123/456";
      const doiUrl = `https://doi.org/${doiPath}`;
      const path = getDOIPath(doiUrl);
      expect(path).toEqual(doiPath);
    });
  });
});
