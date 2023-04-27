import { test } from "@playwright/test";
import { goToWMG } from "../../utils/wmgUtils";
import { isDevStagingProd } from "tests/utils/helpers";
import { SHARED_LINK, SIMPLE_SHARED_LINK } from "tests/common/constants";
import {
  deleteCsvDownloads,
  downloadAndVerifyFiles,
} from "tests/utils/downloadUtils";

const { describe, skip } = test;
describe("SVG download tests", () => {
  //skip(!isDevStagingProd, "WMG BE API does not work locally or in rdev");
  test.only(`Should verify SVG download with single tissue`, async ({
    page,
  }) => {
    await page.pause();
    // set app state
    await goToWMG(page, SIMPLE_SHARED_LINK);

    const tissues = ["blood", "lung"];
    const filterName = "dataset-filter"; // todo: handle multiple filters
    const fileTypes = ["csv"];
    //download and verify svg file
    await downloadAndVerifyFiles(page, filterName, fileTypes, tissues);
  });

  test.afterAll(async () => {
    //delete csv
    //deleteCsvDownloads(`./tests/download`);
  });
});
