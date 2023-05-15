import { test } from "@playwright/test";
import { goToWMG } from "../../utils/wmgUtils";
import {
  SHARED_LINK,
  SIMPLE_SHARED_LINK,
  downLoadPath,
} from "tests/common/constants";
import {
  compareSvg,
  deleteDownloadedFiles,
  downloadAndVerifyFiles,
  subDirectory,
  captureTissueSnapshot,
} from "tests/utils/downloadUtils";
import { isDevStagingProd } from "tests/utils/helpers";

const { describe, skip } = test;

describe("SVG download tests", () => {
  skip(!isDevStagingProd, "WMG BE API does not work locally or in rdev");

  test(`Should verify SVG download without grouping`, async ({ page }) => {
    const tissues = ["blood", "lung"];
    const fileTypes = ["svg"];
    const folder = subDirectory();
    // verify SVG
    for (let i = 0; i < tissues.length; i++) {
      const cellSnapshot = `${downLoadPath}/${folder}/${tissues[i]}.png`;
      const geneSnapshot = `${downLoadPath}/${folder}/gene_${i}.png`;
      // set app state
      await goToWMG(page, SIMPLE_SHARED_LINK);
      //download and verify svg file
      await downloadAndVerifyFiles(page, fileTypes, tissues, folder);
      await captureTissueSnapshot(page, downLoadPath, folder, tissues, i);
      await compareSvg(
        page,
        cellSnapshot,
        geneSnapshot,
        `${folder}/${tissues[i]}.svg`,
        folder,
        tissues[i]
      );
      await deleteDownloadedFiles(`./tests/downloads/${folder}`);
    }
  });

  test(`Should verify SVG download with grouping`, async ({ page }) => {
    const tissues = ["blood"]; // this filter resolves to one tissue
    const fileTypes = ["svg"];
    const folder = subDirectory();
    // verify SVG
    for (let i = 0; i < tissues.length; i++) {
      const cellSnapshot = `${downLoadPath}/${folder}/${tissues[i]}.png`;
      const geneSnapshot = `${downLoadPath}/${folder}/gene_${i}.png`;

      // set app state
      await goToWMG(page, SHARED_LINK);
      //download and verify svg file
      await downloadAndVerifyFiles(page, fileTypes, tissues, folder);
      await captureTissueSnapshot(page, downLoadPath, folder, tissues, i);
      await compareSvg(
        page,
        cellSnapshot,
        geneSnapshot,
        `${folder}/${tissues[i]}.svg`,
        folder,
        tissues[i]
      );
      await deleteDownloadedFiles(`./tests/downloads/${folder}`);
    }
  });
});
