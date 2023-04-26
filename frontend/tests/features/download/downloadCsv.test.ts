import { test } from "@playwright/test";
import {
  goToWMG,
  selectFilterOption,
  selectGroupByOption,
  selectTissueAndGeneOption,
} from "../../utils/wmgUtils";
import { isDevStagingProd } from "tests/utils/helpers";
import { SHARED_LINK } from "tests/common/constants";
import {
  downloadAndVerifyCsv,
  deleteCsvDownloads,
} from "tests/utils/downloadUtils";

const DATASET_FILTER = "dataset-filter";
const { describe, skip } = test;
describe("CSV download tests", () => {
  test.beforeEach(async ({ page }) => {
    // navigate to gene expression page
    await goToWMG(page);
    //select tissue and gene
    await selectTissueAndGeneOption(page);
  });
  //skip(!isDevStagingProd, "WMG BE API does not work locally or in rdev");
  test(`Should verify CSV metadata and header for lung tissue with no filter applied`, async ({
    page,
  }) => {
    //download and verify csv file
    await downloadAndVerifyCsv(page, "no-filter", "lung");
  });
  test(`Should verify CSV metadata and header for blood tissue with no filter applied`, async ({
    page,
  }) => {
    //download and verify csv file
    await downloadAndVerifyCsv(page, "no-filter", "blood");
  });
  test(`Should verify CSV metadata and header for blood tissue with filter applied and group by selected`, async ({
    page,
  }) => {
    //group by
    await selectGroupByOption(page);

    //select a filter
    await selectFilterOption(page, DATASET_FILTER);

    //download and verify csv file
    await downloadAndVerifyCsv(page, DATASET_FILTER, "blood");
  });
  test(`Should verify CSV metadata and header with dataset filter applied and group by selected`, async ({
    page,
  }) => {
  
    // set app state
    await goToWMG(page);

    //select tissue and gene
    await selectTissueAndGeneOption(page);

    //group by
    await selectGroupByOption(page);

    //select a filter
    await selectFilterOption(page, DATASET_FILTER);

    //download and verify csv file
    await downloadAndVerifyCsv(page, DATASET_FILTER, "blood");
  });
  test(`Should verify CSV metadata and header with sex filter applied and group by selected`, async ({
    page,
  }) => {
    // set app state
    await goToWMG(page);

    //select tissue and gene
    await selectTissueAndGeneOption(page);

    //group by
    await selectGroupByOption(page);

    //select a filter
    await selectFilterOption(page, "sex-filter");

    //download and verify csv file
    await downloadAndVerifyCsv(page, "sex-filter", "lung");
  });
  test.afterAll(async () => {
    //delete csv
    deleteCsvDownloads(`./tests/download`);
  });
});
