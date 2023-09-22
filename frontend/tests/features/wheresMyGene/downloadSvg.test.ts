import { subDirectory, verifySvgDownload } from "tests/utils/downloadUtils";
import { test } from "@playwright/test";
import { SHARED_LINK, SIMPLE_SHARED_LINK } from "tests/common/constants";

const { describe, skip } = test;

describe("SVG download tests", () => {
  skip(true, "Skip for now, because image diffing is hard");

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
