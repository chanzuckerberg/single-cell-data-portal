import { subDirectory, verifyImageDownload } from "tests/utils/downloadUtils";
import { test } from "@playwright/test";
import { SHARED_LINK, SIMPLE_SHARED_LINK } from "tests/common/constants";

const { describe } = test;

describe("Image download tests", () => {
  test(`Verify image download without grouping`, async ({ page }) => {
    const tissues = ["blood", "lung"];
    const folder = subDirectory();

    await verifyImageDownload(page, SIMPLE_SHARED_LINK, tissues, folder);
  });

  test(`Verify image download with grouping`, async ({ page }) => {
    const tissues = ["blood"];
    const folder = subDirectory();

    await verifyImageDownload(page, SHARED_LINK, tissues, folder);
  });
});
