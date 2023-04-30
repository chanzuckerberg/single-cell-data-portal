import { test } from "@playwright/test";
import { goToWMG } from "../../utils/wmgUtils";
import { SIMPLE_SHARED_LINK } from "tests/common/constants";
import {
  deleteCsvDownloads,
  downloadAndVerifyFiles,
  subDirectory,
  verifySvg,
} from "tests/utils/downloadUtils";

const { describe } = test;
describe("SVG download tests", () => {
  test.only(`Should verify SVG download with single tissue`, async ({
    page,
  }) => {
    // set app state
    await goToWMG(page, SIMPLE_SHARED_LINK);

    const tissues = ["blood", "lung"];
    const fileTypes = ["svg"];
    const folder = subDirectory();
    //download and verify svg file
    await downloadAndVerifyFiles(page, fileTypes, tissues, folder);

    // verify SVG
    const fixturePath = `./tests/fixtures/svg`;
    verifySvg(tissues, folder, fixturePath);
  });

  // test.afterAll(async () => {
  //   //delete csv
  //   deleteCsvDownloads(`./tests/downloads`);
  // });
});
