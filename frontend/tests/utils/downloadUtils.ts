import { Page, chromium, expect } from "@playwright/test";
import * as fs from "fs";
import readline from "readline";
import AdmZip from "adm-zip";
import { getTestID } from "./selectors";
import { ROUTES } from "src/common/constants/routes";
import { TEST_URL } from "tests/common/constants";
import pixelmatch from "pixelmatch";
import { PNG } from "pngjs";
import sharp from "sharp";
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
const downLoadPath = "./tests/downloads";

export async function verifyCsv(
  page: Page,
  subDirectory: string,
  tissues: string[],
  filterName: string,
  url: string
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

    //verify the number of element in the csv this is the Ui displayed multiplied by the number of genes selected
    expect(csvElementsCount).toEqual(uiElementsCount * 3);

    const options = {
      filterName: filterName,
      data: headers,
    };

    //verify meta data
    await verifyMetadata(page, options, url);

    //verify all the headers are present in the csv
    expect(data[0]).toEqual(expect.arrayContaining(EXPECTED_HEADER));
  });
}

export function subDirectory() {
  return (Math.floor(Math.random() * 90000) + 10000).toString();
}
export async function downloadAndVerifyFiles(
  page: Page,
  fileTypes: string[] = ["png"],
  tissues: string[],
  subDirectory: string
): Promise<void> {
  //download and extract file
  await downloadGeneFile(page, tissues, subDirectory, fileTypes);

  // verify files are available
  fileTypes.forEach((extension: string) => {
    tissues.forEach((tissue) => {
      expect(
        fs.existsSync(`${downLoadPath}/${subDirectory}/${tissue}.${extension}`)
      ).toBeTruthy();
    });
  });
}

export function deleteDownloadedFiles(filePath: string) {
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
      `${downLoadPath}/${subDirectory}/${tissue}.csv`,
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
  options: MetadataVerificationOptions,
  url: string
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
    //this gets the filter displayed on the ui
    text = await getFilterText(page, options.filterName);
  }

  // Check if a match was found and log the result
  if (linkMatch) {
    const link = linkMatch[0];
    //verify link is valid

    expect(link).toEqual(url);

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

  for (const ext of allFileTypes) {
    const checkboxId = `${ext}-checkbox`;
    // uncheck if file type is checked but not
    if (
      (await page.getByTestId(checkboxId).getAttribute("class"))?.includes(
        CHECK
      ) &&
      !fileTypes.includes(ext)
    ) {
      await page.getByTestId(checkboxId).locator("input").click();
    }

    // ensure wanted file types are checked
    if (
      !(await page.getByTestId(checkboxId).getAttribute("class"))?.includes(
        CHECK
      ) &&
      fileTypes.includes(ext)
    ) {
      await page.getByTestId(checkboxId).locator("input").click();
    }
  }
  await page.getByTestId("dialog-download-button").click();
  const download = await downloadPromise;

  // download can be zipped or not depending on number of tissues
  let fileName = `${dirPath}/download.zip`;
  if (fileTypes.length === 1 && tissues.length === 1 && fileTypes[0] == "svg") {
    fileName = `${dirPath}/${tissues[0]}.${fileTypes[0]}`;
  }

  await download.saveAs(fileName);
  //extract zip file
  if (fileName.includes("zip")) {
    const zip = new AdmZip(fileName);
    zip.extractAllTo(dirPath);
  }
}

export const getFilterText = async (page: Page, filterName: string) => {
  const filter_label = `${getTestID(filterName)} [role="button"]`;
  return await page.locator(filter_label).textContent();
};

export function verifyPng(
  dirPath: string,
  tissues: string[],
  subDirectory: string
) {
  tissues.forEach(async (tissue) => {
    // Capture the actual screenshot and compare it with the expected screenshot
    const expectedPng = `./tests/fixtures/download/${subDirectory}/${tissue}.png`;

    const actualPng = `${dirPath}/${tissue}.png`;
    compareImages(expectedPng, actualPng);
  });
}
async function compareImages(imagePath1: string, imagePath2: string) {
  const imageBuffer1 = await fs.promises.readFile(imagePath1);
  const imageBuffer2 = await fs.promises.readFile(imagePath2);

  const image1 = PNG.sync.read(imageBuffer1);
  const image2 = PNG.sync.read(imageBuffer2);

  if (image1.width !== image2.width || image1.height !== image2.height) {
    const maxWidth = Math.max(image1.width, image2.width);
    const maxHeight = Math.max(image1.height, image2.height);

    const resizedImage1 = await sharp(imageBuffer1)
      .resize(maxWidth, maxHeight)
      .png()
      .toBuffer();
    const resizedImage2 = await sharp(imageBuffer2)
      .resize(maxWidth, maxHeight)
      .png()
      .toBuffer();

    const resizedImage1PNG = PNG.sync.read(resizedImage1);
    const resizedImage2PNG = PNG.sync.read(resizedImage2);

    const diffPNG = new PNG({ width: maxWidth, height: maxHeight });
    const numDiffPixels = pixelmatch(
      resizedImage1PNG.data,
      resizedImage2PNG.data,
      diffPNG.data,
      maxWidth,
      maxHeight,
      { threshold: 0.95 }
    );
    expect(numDiffPixels).toBeLessThan(400000);
  }
}
