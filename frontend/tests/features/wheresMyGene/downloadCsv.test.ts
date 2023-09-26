import { test } from "@playwright/test";
import { conditionallyRunTests, goToWMG } from "../../utils/wmgUtils";
import {
  subDirectory,
  downloadAndVerifyFiles,
  verifyCsv,
  deleteDownloadedFiles,
} from "tests/utils/downloadUtils";

import {
  SHARED_LINK_FILTER,
  SHARED_LINK_NO_FILTER,
  SHARED_LINK_NO_GROUP,
} from "tests/common/constants";

const { describe } = test;
describe("CSV download tests", () => {
  conditionallyRunTests({ forceRun: true });

  test(`Should verify CSV metadata and header for lung tissue with no group set`, async ({
    page,
  }) => {
    // set app state
    await goToWMG(page, SHARED_LINK_NO_GROUP);
    const tissues = ["blood", "lung"];
    const fileTypes = ["csv"];
    const folder = subDirectory();
    //download  csv file
    await downloadAndVerifyFiles(page, fileTypes, tissues, folder);
    // verify csv file
    await verifyCsv(page, folder, "disease-filter", SHARED_LINK_NO_GROUP);
    await deleteDownloadedFiles(`./tests/downloads/${folder}`);
  });

  test(`Should verify CSV metadata and header lung and blood tissue with sex filter applied and group by selected`, async ({
    page,
  }) => {
    // set app state
    await goToWMG(page, SHARED_LINK_FILTER);
    const tissues = ["blood", "lung"];
    const fileTypes = ["csv"];
    const folder = subDirectory();
    //download  csv file
    await downloadAndVerifyFiles(page, fileTypes, tissues, folder);

    // verify csv file
    await verifyCsv(page, folder, "sex-filter", SHARED_LINK_FILTER);
    await deleteDownloadedFiles(`./tests/downloads/${folder}`);
  });

  test(`Should verify CSV metadata and header for lung and blood tissue with no filter applied`, async ({
    page,
  }) => {
    // set app state
    await goToWMG(page, SHARED_LINK_NO_FILTER);
    const tissues = ["blood", "lung"];
    const fileTypes = ["csv"];
    const folder = subDirectory();
    //download  csv file
    await downloadAndVerifyFiles(page, fileTypes, tissues, folder);

    // verify csv file
    await verifyCsv(page, folder, "no-filter", SHARED_LINK_NO_FILTER);
    await deleteDownloadedFiles(`./tests/downloads/${folder}`);
  });
});
