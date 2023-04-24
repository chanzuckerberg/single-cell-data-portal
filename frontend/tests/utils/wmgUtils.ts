import { ROUTES } from "src/common/constants/routes";
import { TEST_URL } from "../common/constants";
import { expect, Page } from "@playwright/test";
import { getTestID, getText } from "tests/utils/selectors";
import AdmZip from "adm-zip";
import * as fs from "fs";
import readline from "readline";

export function goToWMG(page: Page) {
  return Promise.all([
    page.waitForResponse(
      (resp: { url: () => string | string[]; status: () => number }) =>
        resp.url().includes("/wmg/v1/filters") && resp.status() === 200
    ),
    page.goto(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`),
  ]);
}

export const selectFilterOption = async (page: Page, filterName: string) => {
  // click the filter at the corner this is done due to the fact that the default click is being intercepted by another element
  await page.getByTestId(filterName).getByRole("button").click();

  // select the first option
  await page.locator("[data-option-index='0']").click();

  // close the pop-up
  await page.getByTestId("dataset-filter").click();

  const filter_label = `${getTestID(filterName)} [role="button"]`;
  // expect the selected filter to be visible
  await expect(page.locator(filter_label)).toBeVisible();

  //wait till loading is complete
  await page.locator(getText("Loading")).waitFor({ state: "hidden" });
};
export const pickOptions = async (page: Page, n: number) => {
  for (let i = 0; i < n; i++) {
    // select the nth option
    await page.locator(`[data-option-index="${i}"]`).click();
  }
};

export const deSelectFilterOption = async (page: Page, filterName: string) => {
  const filter_label = `${getTestID(filterName)} [role="button"]`;
  // expect the selected filter to be visible
  await expect(page.locator(filter_label)).toBeVisible();

  // click the cancel button
  await page.getByTestId("ClearIcon").click();

  // verify the selected filter is not visible
  const visibility = await page.locator(filter_label).isVisible();
  expect(visibility).toBeFalsy();
};

export const selectOption = async (page: Page, filterName: string) => {
  // click the filter
  await page.getByTestId(filterName).click();
};

export const selectTissueAndGeneOption = async (page: Page) => {
  // click Tissue button
  await selectOption(page, "add-tissue-btn");

  //pick the first 2 elements in tissue
  await pickOptions(page, 2);

  // close the pop-up
  await page.keyboard.press("Escape");

  //wait for heatmap to be visible the click action
  await page.locator('[id="heatmap-container-id"]').waitFor();

  // click Gene button
  await selectOption(page, "add-gene-btn");

  //pick the first n elements in tissue
  await pickOptions(page, 3);

  // close the pop-up
  await page.keyboard.press("Escape");

  //wait for gene label to appear
  await page.locator("[data-testid='gene-label-TSPAN6']").waitFor();

  //wait till loading is complete
  await page.locator(getText("Loading")).waitFor({ state: "hidden" });
};

export const checkSourceData = async (page: Page) => {
  //click on source data icon
  await page.locator('[data-testid="source-data-button"]').click();

  // number of elemet displayed on source data
  const n = await page.locator('[data-testid="source-data-list"] a').count();
  // close the pop-up
  await page.keyboard.press("Escape");
  return n;
};
export const checkPlotSize = async (page: Page) => {
  //get the number of rows on the data plot
  const n = await page.locator('[data-zr-dom-id*="zr"]').count();
  let sumOfHeights = 0;
  for (let i = 0; i < n; i++) {
    const row = await page.locator('[data-zr-dom-id*="zr"]').nth(i);

    const height = await row.getAttribute("height");

    if (height !== null) {
      sumOfHeights += parseInt(height);
    }
  }
  return sumOfHeights;
};
export const downloadCsv = async (page: Page, fileFactor: string) => {
  const zipFilePath = "./tests/utils/download.zip";
  const extractDirPath = `./tests/utils/download/${fileFactor}`;
  const CHECK = "Mui-checked";
  //wait for download file
  const downloadPromise = page.waitForEvent("download");

  //click the download icon

  await page.getByTestId("download-button").click();
  const checkboxClassPng = await page
    .getByTestId("png-checkbox")
    .getAttribute("class");

  if (checkboxClassPng && checkboxClassPng.includes(CHECK)) {
    await page.getByTestId("png-checkbox").click();
  }
  const checkboxClassSvg = await page
    .getByTestId("svg-checkbox")
    .getAttribute("class");

  if (checkboxClassSvg && checkboxClassSvg.includes(CHECK)) {
    await page.getByTestId("svg-checkbox").click();
  }

  const checkboxClassCsv = await page
    .getByTestId("csv-checkbox")
    .getAttribute("class");

  if (checkboxClassCsv && !checkboxClassCsv.includes(CHECK)) {
    await page.getByTestId("csv-checkbox").click();
  }

  await page.getByTestId("dialog-download-button").click();
  const download = await downloadPromise;

  // Save downloaded file in a directory
  await download.saveAs(zipFilePath);

  //extract zip file
  const zip = new AdmZip(zipFilePath);
  zip.extractAllTo(extractDirPath);
};

export const getFilterText = async (page: Page, filterName: string) => {
  const filter_label = `${getTestID(filterName)} [role="button"]`;
  return await page.locator(filter_label).textContent();
};
export const verifyMetadata = async (
  page: Page,
  filterName: string,
  data: string[],
  noSelectionText = "No selection"
) => {
  //verify the date is valid
  const dateString = data[0].substring(14);
  const date = new Date(dateString);

  // Check if the resulting date is valid
  expect(!isNaN(date.getTime())).toBe(true);

  expect(
    data[1].includes(
      "# We regularly expand our single cell data corpus to improve results. Downloaded data and figures may differ in the future."
    )
  ).toBeTruthy();

  //check if the 3 column contains the gene expression link
  expect(
    data[2].includes("https://localhost:3000/gene-expression")
  ).toBeTruthy();

  // Extract the link using a regular expression
  const linkRegex = /https?:\/\/[^\s]+/;
  const linkMatch = data[2].match(linkRegex);
  let text: string | null = "";
  if (filterName !== "no-filter") {
    text = await getFilterText(page, filterName);
  }

  // Check if a match was found and log the result
  if (linkMatch) {
    const link = linkMatch[0];
    //verify link is valid

    await page.goto(link);

    // wait until the new page fully loads
    await page.waitForLoadState();
    // expect the header on the new page to be visible
    expect(page.getByTestId("download-button")).toBeVisible();
  }

  //verify the text displayed in the filter section of meta data
  switch (filterName) {
    case "dataset-filter":
      expect(
        data[3].includes("# Dataset Filter Values: " + `${text}`)
      ).toBeTruthy();
      expect(
        data[4].includes(`# Disease Filter Values: ${noSelectionText}`)
      ).toBeTruthy();
      expect(
        data[5].includes(
          `# Self-Reported Ethnicity Filter Values: ${noSelectionText}`
        )
      ).toBeTruthy();
      expect(
        data[6].includes(`# Sex Filter Values: ${noSelectionText}`)
      ).toBeTruthy();
      break;
    case "disease-filter":
      expect(data[4].includes(`# Disease Filter Values: ${text}`)).toBeTruthy();
      expect(
        data[3].includes(`# Dataset Filter Values: ${noSelectionText}`)
      ).toBeTruthy();
      expect(
        data[5].includes(
          `# Self-Reported Ethnicity Filter Values: ${noSelectionText}`
        )
      ).toBeTruthy();
      expect(
        data[6].includes(`# Sex Filter Values: ${noSelectionText}`)
      ).toBeTruthy();
      break;
    case "self-reported-ethnicity-filter":
      expect(
        data[5].includes(`# Self-Reported Ethnicity Filter Values: ${text}`)
      ).toBeTruthy();
      expect(
        data[3].includes(`# Dataset Filter Values: ${noSelectionText}`)
      ).toBeTruthy();
      expect(
        data[4].includes(`# Disease Filter Values: ${noSelectionText}`)
      ).toBeTruthy();
      expect(
        data[6].includes(`# Sex Filter Values: ${noSelectionText}`)
      ).toBeTruthy();
      break;
    case "sex-filter":
      expect(data[6].includes(`# Sex Filter Values: ${text}`)).toBeTruthy();
      expect(
        data[3].includes(`# Dataset Filter Values: ${noSelectionText}`)
      ).toBeTruthy();
      expect(
        data[4].includes(`# Disease Filter Values: ${noSelectionText}`)
      ).toBeTruthy();
      expect(
        data[5].includes(
          `# Self-Reported Ethnicity Filter Values: ${noSelectionText}`
        )
      ).toBeTruthy();
      break;
    case "no-filter":
      expect(
        data[6].includes(`# Sex Filter Values: ${noSelectionText}`)
      ).toBeTruthy();
      expect(
        data[3].includes(`# Dataset Filter Values: ${noSelectionText}`)
      ).toBeTruthy();
      expect(
        data[4].includes(`# Disease Filter Values: ${noSelectionText}`)
      ).toBeTruthy();
      expect(
        data[5].includes(
          `# Self-Reported Ethnicity Filter Values: ${noSelectionText}`
        )
      ).toBeTruthy();
      break;
    default:
      throw new Error(`Invalid filter name: ${filterName}`);
  }
};

export const deleteCsvFile = (filePath: string): void => {
  fs.rmdir(filePath, { recursive: true }, (err) => {
    if (err) {
      console.error(`Error deleting folder: ${err}`);
    } else {
      console.log("Folder deleted successfully");
    }
  });
};

export interface CsvMetadata {
  rowCount: number;
  data: string[][];
  headers: string[];
}

export const getCsvMetadata = (
  tissue: string,
  fileFactor: string
): Promise<CsvMetadata> => {
  return new Promise((resolve, reject) => {
    // Open the CSV file for reading
    const fileStream = fs.createReadStream(
      `./tests/utils/download/${fileFactor}/${tissue}.csv`,
      { encoding: "utf8" }
    );

    // Create a readline interface for the file stream
    const rl = readline.createInterface({ input: fileStream });

    // Create counters and arrays to store the parsed data
    let rowCount = -1;
    const data: string[][] = [];
    const headers: string[] = [];

    // Listen for 'line' events emitted by the readline interface
    rl.on("line", (line) => {
      const row = line.split(",");
      if (row.length > 1) {
        data.push(row);
        rowCount++;
      } else if (row.length === 1) {
        headers.push(row[0]);
      }
    });

    // Listen for the 'close' event to know when the parsing is complete
    rl.on("close", () => {
      const csvMetadata: CsvMetadata = { rowCount, data, headers };
      resolve(csvMetadata);
    });

    // Listen for any errors during reading and parsing the file
    rl.on("error", (err) => {
      reject(err);
    });
  });
};
