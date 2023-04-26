import { test } from "@playwright/test";
import { goToWMG } from "../../utils/wmgUtils";
import { isDevStagingProd } from "tests/utils/helpers";
import { SHARED_LINK } from "tests/common/constants";
import {
  deleteCsvDownloads,
  downloadAndVerifyCsv,
} from "tests/utils/downloadUtils";

const { describe, skip } = test;
describe("SVG download tests", () => {
  skip(!isDevStagingProd, "WMG BE API does not work locally or in rdev");
  test(`Should verify SVG download with filters applied`, async ({ page }) => {
    // set app state
    await goToWMG(page, SHARED_LINK);

    //download and verify csv file
    await downloadAndVerifyCsv(page, "dataset-filter");
  });

  test.afterAll(async () => {
    //delete csv
    deleteCsvDownloads(`./tests/download`);
  });
});
