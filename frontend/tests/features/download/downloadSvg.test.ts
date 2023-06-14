import { subDirectory, verifySvgDownload } from "tests/utils/downloadUtils";
import { test } from "@playwright/test";
import { isDevStagingProd } from "tests/utils/helpers";
import { SHARED_LINK, SIMPLE_SHARED_LINK } from "tests/common/constants";

const { describe, skip } = test;

describe("SVG download tests", () => {
  skip(!isDevStagingProd, "WMG BE API does not work locally or in rdev");

  test(`Should verify SVG download without grouping`, async ({ page }) => {
    const tissues = ["blood", "lung"];
    const folder = subDirectory();

    await verifySvgDownload(page, SIMPLE_SHARED_LINK, tissues, folder);
  });

  test(`Should verify SVG download with grouping`, async ({ page }) => {
    const tissues = ["blood"];
    const folder = subDirectory();

    await verifySvgDownload(page, SHARED_LINK, tissues, folder);
  });
});
