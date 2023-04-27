import { Page, expect } from "@playwright/test";
import * as fs from "fs";
import readline from "readline";
import AdmZip from "adm-zip";
import { getTestID } from "./selectors";
import { ROUTES } from "src/common/constants/routes";
import { TEST_URL } from "tests/common/constants";

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
const downLoadPath = "./tests/download/";

export async function verifyCsv(
  page: Page,
  subDirectory: string,
  tissues: string[],
  filterName: string
): Promise<void> {
  tissues.forEach(async (tissue) => {
    const metadata = await getCsvMetadata(tissue, subDirectory);

    // extract the headers and data arrays from the metadata object
    // put all the headers in an array
    const headers = metadata.headers;
    const data = metadata.data;

    //get number of element in csv
    const csvElementsCount = metadata.rowCount;

    //get number of element displayed in ui
    const uiElementsCount = await page
      .locator(
        `[data-testid="cell-type-labels-${tissue}"] [data-testid="cell-type-label-count"]`
      )
      .count();

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
  });
}

export async function downloadAndVerifyFiles(
  page: Page,
  filterName: string,
  fileTypes: string[] = ["png"],
  tissues: string[]
): Promise<void> {
  // generate sub folder
  const subDirectory: string = (
    Math.floor(Math.random() * 90000) + 10000
  ).toString();

  //download and extract file
  await downloadGeneFile(page, tissues, subDirectory, fileTypes);

  // verify files are available
  // fileTypes.forEach((extension: string) => {
  //   expect(
  //     fs.existsSync(`${downLoadPath}/${subDirectory}/${extension}`)
  //   ).toBeTruthy();
  // });

  // verify CSV
  if (fileTypes.includes("csv")) {
    //verify meta data
    await verifyCsv(page, subDirectory, tissues, filterName);
  }
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
  subDirectory: string
): Promise<CsvMetadata> => {
  return new Promise((resolve, reject) => {
    // Open the CSV file for reading
    const fileStream = fs.createReadStream(
      `${downLoadPath}${subDirectory}/${tissue}.csv`,
      { encoding: "utf8" }
    );

    // Create a readline interface for the file stream
    const csvFileInterface = readline.createInterface({ input: fileStream });

    // Create counters and arrays to store the parsed data
    let rowCount = -1;
    const data: string[][] = [];
    const headers: string[] = [];

    // Listen for 'line' events emitted by the readline interface
    csvFileInterface.on("line", (line) => {
      const row = line.split(",");
      if (row.length > 1) {
        data.push(row);
        rowCount++;
      } else if (row.length === 1) {
        headers.push(row[0]);
      }
    });

    // Listen for the 'close' event to know when the parsing is complete
    csvFileInterface.on("close", () => {
      const csvMetadata: CsvMetadata = { rowCount, data, headers };
      resolve(csvMetadata);
    });

    // Listen for any errors during reading and parsing the file
    csvFileInterface.on("error", (err) => {
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
    options.data[2].includes(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`)
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

export async function downloadGeneFile(
  page: Page,
  tissues: string[],
  subDirectory: string,
  fileTypes: string[] = ["png"]
): Promise<void> {
  const allFileTypes = ["csv", "png", "svg"];
  const dirPath = `${downLoadPath}/${subDirectory}`;

  //wait for download file
  const downloadPromise = page.waitForEvent("download");

  //click the download icon
  await page.getByTestId("download-button").click();
  const CHECK = "Mui-checked";
  // uncheck and check file types as needed
  for (const _fileType in allFileTypes) {
    const checkboxId = `${allFileTypes[_fileType]}-checkbox`;
    // uncheck unwanted file type
    if (
      (await page.getByTestId(checkboxId).getAttribute("class"))?.includes(
        CHECK
      ) &&
      !fileTypes.includes(allFileTypes[_fileType])
    ) {
      await page.getByTestId(checkboxId).click();
    }

    // ensure wanted file types are checked
    if (
      !(await page.getByTestId(checkboxId).getAttribute("class"))?.includes(
        CHECK
      ) &&
      fileTypes.includes(allFileTypes[_fileType])
    ) {
      await page.getByTestId(checkboxId).click();
    }
  }

  await page.getByTestId("dialog-download-button").click();
  const download = await downloadPromise;

  // download can be zipped or not depending on number of tissues
  // if (fileTypes.length === 1 || tissues.length === 1) {
  //   const fileName = `${dirPath}/download.${allFileTypes[0]}`;
  //   await download.saveAs(fileName);
  // } else {
  const fileName = `${dirPath}/download.zip`;
  await download.saveAs(fileName);
  //extract zip file
  if (fileName.includes("zip")) {
    const zip = new AdmZip(fileName);
    zip.extractAllTo(dirPath);
    // }
  }
}

export const getFilterText = async (page: Page, filterName: string) => {
  const filter_label = `${getTestID(filterName)} [role="button"]`;
  return await page.locator(filter_label).textContent();
};
