import { expect, test } from "@playwright/test";
import { goToWMG, selectTissueAndGeneOption } from "../../utils/wmgUtils";
import { isDevStagingProd } from "tests/utils/helpers";
import {
  deleteDownloadedFiles,
  downloadAndVerifyFiles,
  subDirectory,
  verifyPng,
} from "tests/utils/downloadUtils";
import pixelmatch from "pixelmatch";
import fs from "fs";
import { PNG } from "pngjs";
import {
  SHARED_LINK_NO_FILTER,
  SHARED_LINK_NO_GROUP,
} from "tests/common/constants";
const downLoadPath = "./tests/downloads";
const { describe, skip } = test;
describe("CSV download tests", () => {
  //skip(!isDevStagingProd, "WMG BE API does not work locally or in rdev");

  test.only(`Should verify png for lung and blood tissue with no group set`, async ({
    page,
  }) => {
    // set app state
    await goToWMG(page, SHARED_LINK_NO_FILTER);
    const tissues = ["blood"];
    const fileTypes = ["png"];
    const folder = subDirectory();
    const dirPath = `${downLoadPath}/${folder}`;
    //download  csv file
    await downloadAndVerifyFiles(page, fileTypes, tissues, folder);

    verifyPng(dirPath, tissues, "no-filter");
  });

  test.afterAll(async () => {
    //delete csv
   // deleteDownloadedFiles(`./tests/downloads`);
  });
});
