import { test } from "@playwright/test";
import { goToWMG, selectTissueAndGeneOption } from "../../utils/wmgUtils";
import * as fs from "fs";
import AdmZip from "adm-zip";
import { isDevStagingProd } from "tests/utils/helpers";
import readline from "readline";

const { describe, skip } = test;

describe("Left side bar", () => {
  // skip(!isDevStagingProd, "WMG BE API does not work locally or in rdev");
  test.beforeEach(async ({ page }) => {
    // navigate to gene expression page
    await goToWMG(page);

    //select tissue and gene
    await selectTissueAndGeneOption(page);
  });

  test.only("check download", async ({ page }) => {
    // await page.pause();
    const downloadPromise = page.waitForEvent("download");
    await page.getByTestId("download-button").click();
    await page.getByTestId("png-checkbox").click();
    await page.getByTestId("csv-checkbox").click();
    await page.getByTestId("dialog-download-button").click();
    const download = await downloadPromise;
    // Wait for the download process to complete
    console.log(await download.path());
    // Save downloaded file somewhere
    await download.saveAs("./tests/utils/download.zip");
    const zipFilePath = "./tests/utils/download.zip";
    const extractDirPath = "./tests/utils/download";
    const zip = new AdmZip(zipFilePath);
    zip.extractAllTo(extractDirPath);
    const csvFilePath = "./tests/utils/download/blood.csv";

    // Open the CSV file for reading
    const fileStream = fs.createReadStream("./tests/utils/download/blood.csv", {
      encoding: "utf8",
    });

    // Create a readline interface for the file stream
    const rl = readline.createInterface({ input: fileStream });

    // Create an empty array to store the parsed data
    let data: string[]=[];

    // Listen for 'line' events emitted by the readline interface
    rl.on("line", (line) => {
      const row = line.split(",");
      if (row.length == 1) {
        data.push(row[0]);
      }
    });

    // Listen for the 'close' event to know when the parsing is complete
    rl.on("close", () => {
      const outputStream = fs.createWriteStream(
        "./tests/utils/download/file.json"
      );
      outputStream.write(JSON.stringify(data));
      outputStream.end();
    });
  });
});
