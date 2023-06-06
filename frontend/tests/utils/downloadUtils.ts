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

export async function deleteDownloadedFiles(filePath: string) {
  fs.rmdir(filePath, { recursive: true }, (err) => {
    if (err) {
      console.error(`Error deleting folder: ${err}`);
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
interface MetadataVerificationOptions {
  filterName: string;
  data: string[];
  noSelectionText?: string;
}

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
  const fileName = `${dirPath}/download.zip`;
  // if (
  //   fileTypes.length === 1 &&
  //   tissues.length === 1 &&
  //   fileTypes[0] !== "png"
  // ) {
  //   fileName = `${dirPath}/${tissues[0]}.${fileTypes[0]}`;
  // }

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

export function verifyPng(page: Page, dirPath: string, tissues: string[]) {
  tissues.forEach(async () => {
    // Specify the input image file path
    const inputImagePath = `${dirPath}/blood.png`;

    // Specify the output image file path
    const outputImagePath = `${dirPath}/blood_cropped.png`;
    // Specify the height of the region to be removed from the top
    const startHeight = 567;
    // Perform the extraction operation
    sharp(inputImagePath)
      .metadata()
      .then((metadata) => {
        const imageHeight = metadata.height || 0;
        const extractHeight = imageHeight - startHeight;

        const extractRegion = {
          top: startHeight,
          left: 0,
          width: metadata.width || 0,
          height: extractHeight,
        };

        sharp(inputImagePath)
          .extract(extractRegion)
          .toFile(outputImagePath)
          .then(() => {
            console.log("Image processed successfully");
          })
          .catch((err) => {
            console.error("Error while processing image:", err);
          });
      })
      .catch((err) => {
        console.error("Error while retrieving image metadata:", err);
      });

    const actualPng = outputImagePath;
    // Capture the actual screenshot and compare it with the expected screenshot
    const expectedPng = `${dirPath}/blood_actual.png`;
    await page.waitForTimeout(1000);
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
    console.log(numDiffPixels);
    expect(numDiffPixels).toBeLessThan(400000);
  }
}
export async function compareSvg(
  webCellImage: string,
  webGeneImage: string,
  svgFile: string
): Promise<void> {
  const browser = await chromium.launch();
  const browserContext = await browser.newContext();
  const page = await browserContext.newPage();
  await page.goto(svgFile);
  expect(page.locator("svg").locator("svg").nth(3)).toMatchSnapshot(
    webCellImage
  );
  expect(page.locator("svg").locator("svg").nth(4)).toMatchSnapshot(
    webCellImage
  );
  expect(await page.screenshot()).toMatchSnapshot(webGeneImage);
  await browser.close();
}
export async function getActualImage(
  page: Page,
  dirPath: string
): Promise<void> {
  const geneCanvasId = '[data-zr-dom-id*="zr"]';
  //close popup
  await page
    .locator('[data-testid="newsletter-modal-banner-wrapper"] button')
    .click();
  // Get the start element
  const startElement = await page.$("#heatmap-container-id");
  if (!startElement) {
    console.error("Start element not found");
    await page.close();
    return;
  }

  // Get the stop element
  const stopElement = await page.locator(geneCanvasId).nth(0);
  if (!stopElement) {
    console.error("Stop element not found");
    await page.close();
    return;
  }

  // Get the bounding boxes of the elements
  const startBoundingBox = await startElement.boundingBox();
  const stopBoundingBox = await stopElement.boundingBox();

  if (!startBoundingBox || !stopBoundingBox) {
    console.error("Bounding box not found");
    await page.close();
    return;
  }

  // Calculate the starting and stopping coordinates
  const startX = startBoundingBox.x;
  const startY = startBoundingBox.y;
  const stopX = stopBoundingBox.x + stopBoundingBox.width;
  const stopY = stopBoundingBox.y + stopBoundingBox.height;

  // Calculate the width and height of the desired screenshot area
  const width = stopX - startX;
  const height = stopY - startY;
  // Calculate the viewport dimensions to accommodate the elements
  const viewportWidth = Math.max(
    startBoundingBox.x + startBoundingBox.width,
    stopBoundingBox.x + stopBoundingBox.width
  );
  const viewportHeight = Math.max(
    startBoundingBox.y + startBoundingBox.height,
    stopBoundingBox.y + stopBoundingBox.height
  );

  // Set the viewport size
  await page.setViewportSize({
    width: viewportWidth,
    height: viewportHeight,
  });

  // Take a screenshot of the desired area
  await page.screenshot({
    path: `${dirPath}/blood_actual.png`,
    clip: { x: startX, y: startY, width, height },
  });
}
