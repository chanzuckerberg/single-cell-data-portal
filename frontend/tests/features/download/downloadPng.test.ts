import { test } from "@playwright/test";
import { goToWMG } from "../../utils/wmgUtils";
import { isDevStagingProd } from "tests/utils/helpers";
import {
  deleteDownloadedFiles,
  downloadAndVerifyFiles,
  subDirectory,
  verifyPng,
} from "tests/utils/downloadUtils";
import {
  SHARED_LINK_FILTER,
  SHARED_LINK_NO_FILTER,
  SHARED_LINK_NO_GROUP,
} from "tests/common/constants";
const downLoadPath = "./tests/downloads";
const { describe, skip } = test;
describe("PNG download tests", () => {
  skip(!isDevStagingProd, "WMG BE API does not work locally or in rdev");

  test(`Should verify png for lung and blood tissue with no filter set`, async ({
    page,
  }) => {
    // set app state
    await goToWMG(page, SHARED_LINK_NO_FILTER);
    const tissues = ["blood", "lung"];
    const fileTypes = ["png"];
    const folder = subDirectory();
    const dirPath = `${downLoadPath}/${folder}`;
    //download file
    await downloadAndVerifyFiles(page, fileTypes, tissues, folder);
    verifyPng(dirPath, tissues, "no-filter");
    deleteDownloadedFiles(dirPath);
  });

  test(`Should verify png data for lung and blood tissue with sex filter applied and group by selected`, async ({
    page,
  }) => {
    // set app state
    await goToWMG(page, SHARED_LINK_FILTER);
    const tissues = ["blood", "lung"];
    const fileTypes = ["png"];
    const folder = subDirectory();
    const dirPath = `${downLoadPath}/${folder}`;
    //download file
    await downloadAndVerifyFiles(page, fileTypes, tissues, folder);
    verifyPng(dirPath, tissues, "filter");
    deleteDownloadedFiles(dirPath);
  });
  test(`Should verify png for lung and blood tissue with no group set`, async ({
    page,
  }) => {
    // set app state
    await goToWMG(page, SHARED_LINK_NO_GROUP);
    const tissues = ["blood", "lung"];
    const fileTypes = ["png"];
    const folder = subDirectory();
    const dirPath = `${downLoadPath}/${folder}`;
    //download file
    await downloadAndVerifyFiles(page, fileTypes, tissues, folder);
    verifyPng(dirPath, tissues, "no-group");
    deleteDownloadedFiles(dirPath);
  });
});
