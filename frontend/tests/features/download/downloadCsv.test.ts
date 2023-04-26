import { test } from "@playwright/test";
import {
  deleteCsvDownloads,
  goToWMG,
} from "../../utils/wmgUtils";
import { isDevStagingProd } from "tests/utils/helpers";
import { SHARED_LINK } from "tests/common/constants";
import { downloadAndVerifyCsv } from "tests/utils/downloadUtils";

const filters = [
  ["dataset-filter"],
  ["disease-filter"],
  ["self-reported-ethnicity-filter"],
  ["sex-filter"],
  ["no-filter"],
];

const { describe, skip } = test;
describe("CSV download tests", () => {
  skip(!isDevStagingProd, "WMG BE API does not work locally or in rdev");
  test(`Should verify CSV metadata and header with filters applied`, async ({
    page,
  }) => {
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
