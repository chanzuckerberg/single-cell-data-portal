import { Page, expect } from "@playwright/test";
import * as fs from "fs";
import readline from "readline";
import AdmZip from "adm-zip";
import { getTestID } from "./selectors";

const BLOOD_TISSUE_COUNT =
  '[data-testid="cell-type-labels-blood"] [data-testid="cell-type-label-count"]';

const EXPECTED_HEADER = [
  "Tissue",
  "Cell Type",
  "Cell Count",
  "Tissue Composition",
  "Gene Symbol",
  "Expression",
  '"Expression',
  ' Scaled"',
  "Number of Cells Expressing Genes",
];
export async function downloadAndVerifyCsv(
  page: Page,
  filterName: string
): Promise<void> {
  // generate sub folder
  const randomNumber: number = Math.floor(Math.random() * 90000) + 10000;
  const subDirectory: string = randomNumber.toString();

  //download and extract the csv file
  await downloadCsv(page, subDirectory);
  const metadata = await getCsvMetadata("blood", subDirectory);

  // extract the headers and data arrays from the metadata object
  // put all the headers in an array
  const headers = metadata.headers;
  const data = metadata.data;

  //get number of element in csv
  const csvElementsCount = metadata.rowCount;

  //get number of element displayed in ui
  const uiElementsCount = await page.locator(BLOOD_TISSUE_COUNT).count();

  //verify the number of element in the csv
  expect(csvElementsCount).toEqual(uiElementsCount * 3);

  const options = {
    filterName: filterName,
    data: headers,
  };

  //verify meta data
  await verifyMetadata(page, options);

  //verify all the headers are present in the csv
  expect(data[0]).toEqual(expect.arrayContaining(EXPECTED_HEADER));
}

export function deleteCsvDownloads(filePath: string) {
  fs.rmdir(filePath, { recursive: true }, (err) => {
    if (err) {
      console.error(`Error deleting folder: ${err}`);
    } else {
      console.log("Folder deleted successfully");
    }
  });
}

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
      `./tests/download/${fileFactor}/${tissue}.csv`,
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

interface MetadataVerificationOptions {
  filterName: string;
  data: string[];
  noSelectionText?: string;
}

export const verifyMetadata = async (
  page: Page,
  options: MetadataVerificationOptions
) => {
  //verify the date is valid
  const dateString = options.data[0].substring(14);
  const date = new Date(dateString);

  // Check if the resulting date is valid
  expect(!isNaN(date.getTime())).toBe(true);

  expect(
    options.data[1].includes(
      "# We regularly expand our single cell data corpus to improve results. Downloaded data and figures may differ in the future."
    )
  ).toBeTruthy();

  //check if the 3 column contains the gene expression link
  expect(
    options.data[2].includes("https://localhost:3000/gene-expression")
  ).toBeTruthy();

  // Extract the link using a regular expression
  const linkRegex = /https?:\/\/[^\s]+/;
  const linkMatch = options.data[2].match(linkRegex);
  let text: string | null = "";
  if (options.filterName !== "no-filter") {
    text = await getFilterText(page, options.filterName);
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

  verifyFilterValues(
    options.data,
    options.filterName,
    text,
    options.noSelectionText
  );
};

const verifyFilterValues = (
  data: string[],
  filterValue: string,
  filterText: string | null,
  noSelectionText = "No selection"
) => {
  let datasetFilterValue = noSelectionText;
  let diseaseFilterValue = noSelectionText;
  let ethnicityFilterValue = noSelectionText;
  let sexFilterValue = noSelectionText;

  switch (filterValue) {
    case "dataset-filter":
      datasetFilterValue = filterText || noSelectionText;
      break;
    case "disease-filter":
      diseaseFilterValue = filterText || noSelectionText;
      break;
    case "self-reported-ethnicity-filter":
      ethnicityFilterValue = filterText || noSelectionText;
      break;
    case "sex-filter":
      sexFilterValue = filterText || noSelectionText;
      break;
    case "no-filter":
      break;
    default:
      throw new Error(`Invalid filter name: ${filterValue}`);
  }

  expect(
    data[3].includes(`# Dataset Filter Values: ${datasetFilterValue}`)
  ).toBeTruthy();
  expect(
    data[4].includes(`# Disease Filter Values: ${diseaseFilterValue}`)
  ).toBeTruthy();
  expect(
    data[5].includes(
      `# Self-Reported Ethnicity Filter Values: ${ethnicityFilterValue}`
    )
  ).toBeTruthy();
  expect(
    data[6].includes(`# Sex Filter Values: ${sexFilterValue}`)
  ).toBeTruthy();
};

export const downloadCsv = async (page: Page, fileFactor: string) => {
  const zipFilePath = "./tests/download.zip";
  const extractDirPath = `./tests/download/${fileFactor}`;
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
