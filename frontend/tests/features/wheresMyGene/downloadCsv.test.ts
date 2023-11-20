import { expect } from "@playwright/test";
import { goToWMG } from "../../utils/wmgUtils";
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
import { test } from "tests/common/test";

const { describe } = test;
describe("CSV download tests", () => {
  test(`Should verify CSV metadata and header for lung tissue with no group set`, async ({
    page,
  }) => {
    // set app state
    await goToWMG(page, SHARED_LINK_NO_GROUP);
    await expect(page.locator("canvas")).not.toHaveCount(0, { timeout: 20000 });

    const tissues = ["blood", "lung"];
    const fileTypes = ["csv"];
    const folder = subDirectory();

    // download  csv file
    await downloadAndVerifyFiles(page, fileTypes, tissues, folder);
    // verify csv file
    await verifyCsv({
      page,
      tissues,
      subDirectory: folder,
      filterName: "disease-filter",
      url: SHARED_LINK_NO_GROUP,
    });
    await deleteDownloadedFiles(`./tests/downloads/${folder}`);
  });

  test(`Should verify CSV metadata and header lung and blood tissue with sex filter applied and group by selected`, async ({
    page,
  }) => {
    // set app state
    await goToWMG(page, SHARED_LINK_FILTER);
    await expect(page.locator("canvas")).not.toHaveCount(0, { timeout: 20000 });

    const tissues = ["blood", "lung"];
    const fileTypes = ["csv"];
    const folder = subDirectory();
    //download  csv file
    await downloadAndVerifyFiles(page, fileTypes, tissues, folder);

    // verify csv file
    await verifyCsv({
      page,
      tissues,
      subDirectory: folder,
      filterName: "sex-filter",
      url: SHARED_LINK_FILTER,
    });
    await deleteDownloadedFiles(`./tests/downloads/${folder}`);
  });

  test(`Should verify CSV metadata and header for lung and blood tissue with no filter applied`, async ({
    page,
  }) => {
    // set app state
    await goToWMG(page, SHARED_LINK_NO_FILTER);
    await expect(page.locator("canvas")).not.toHaveCount(0, { timeout: 20000 });

    const tissues = ["blood", "lung"];
    const fileTypes = ["csv"];
    const folder = subDirectory();
    //download  csv file
    await downloadAndVerifyFiles(page, fileTypes, tissues, folder);

    // verify csv file
    await verifyCsv({
      page,
      tissues,
      subDirectory: folder,
      filterName: "no-filter",
      url: SHARED_LINK_NO_FILTER,
    });
    await deleteDownloadedFiles(`./tests/downloads/${folder}`);
  });
});
