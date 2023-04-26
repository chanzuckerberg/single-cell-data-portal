/* eslint-disable @typescript-eslint/no-explicit-any */
import { Page, expect, test } from "@playwright/test";
import { LATEST_SHARE_LINK_VERSION } from "src/views/WheresMyGene/components/GeneSearchBar/components/ShareButton/utils";
import { TEST_URL } from "tests/common/constants";
import { goToWMG, addTissuesAndGenes } from "tests/utils/wmgUtils";

const { describe } = test;
const SHARE_BUTTON = "share-button";
const tissues = ["blood", "lung"];
const genes = ["DPM1", "TNMD", "TSPAN6"];
const dataSets = [
  {
    id: "c874f155-9bf9-4928-b821-f52c876b3e48",
    text: "49 years old male - Fresh PBMCs (1 day post-intubation)",
  },
  {
    id: "db59611b-42de-4035-93aa-1ed39f38b467",
    text: "49 years old male - Fresh PBMCs (2 days post-intubation)",
  },
  {
    id: "eeacb0c1-2217-4cf6-b8ce-1f0fedf1b569",
    text: "49 years old male - Fresh PBMCs (3 days post-intubation)",
  },
  {
    id: "881fe679-c6e0-45a3-9427-c4e81be6921f",
    text: "66 years old female - Fresh PBMCs (2 days post-intubation)",
  },
  {
    id: "ea786a06-5855-48b7-80d7-0313a21a2044",
    text: "66 years old female - Fresh PBMCs (3 days post-intubation)",
  },
  {
    id: "456e8b9b-f872-488b-871d-94534090a865",
    text: "Single-cell atlas of peripheral immune response to SARS-CoV-2 infection",
  },
];
const diseases = [{ id: "MONDO:0100096", text: "COVID-19" }];
const ethnicities = ["unknown"];
const sexes = [
  { id: "PATO:0000383", text: "female" },
  { id: "PATO:0000384", text: "male" },
];
const COMPARE = "disease";
const initialState =
  "https://localhost:3000/gene-expression?datasets=c874f155-9bf9-4928-b821-f52c876b3e48%2Cdb59611b-42de-4035-93aa-1ed39f38b467%2Ceeacb0c1-2217-4cf6-b8ce-1f0fedf1b569%2C881fe679-c6e0-45a3-9427-c4e81be6921f%2Cea786a06-5855-48b7-80d7-0313a21a2044%2C456e8b9b-f872-488b-871d-94534090a865&diseases=MONDO%3A0100096&ethnicities=unknown&sexes=PATO%3A0000383%2CPATO%3A0000384&tissues=blood%2Clung&genes=DPM1%2CTNMD%2CTSPAN6&ver=2&compare=disease";
describe("Share link tests", () => {
  test("Should share link with single tissue and single gene", async ({
    page,
    browserName,
  }) => {
    skipFirefox(browserName);
    const _tissues = ["blood"];
    const _genes = ["SCYL3"];

    //set up sate
    await setupStateAndCopyShareLink(page, _tissues, _genes);
    //verify link
    await verifyShareLink(page, LATEST_SHARE_LINK_VERSION, _tissues, _genes);
  });

  test("Should share link with multiple tissues and multiple genes", async ({
    page,
    browserName,
  }) => {
    skipFirefox(browserName);

    await setupStateAndCopyShareLink(page, tissues, genes);

    // verify link
    await verifyShareLink(page, LATEST_SHARE_LINK_VERSION, tissues, genes);
  });

  test("Should generate share link with correct format for all query param types", async ({
    page,
    browserName,
  }) => {
    skipFirefox(browserName);

    // prepare state
    await goToWMG(page, initialState);

    // verify link parameters and app state
    await verifyShareLink(
      page,
      LATEST_SHARE_LINK_VERSION,
      tissues,
      genes,
      dataSets,
      sexes,
      diseases,
      ethnicities,
      COMPARE
    );
  });
});

async function setupStateAndCopyShareLink(
  page: Page,
  tissues: string[],
  genes: string[]
) {
  await goToWMG(page);
  await expect(page.getByTestId(SHARE_BUTTON)).toBeDisabled();

  // add tissues and genes
  await addTissuesAndGenes(page, tissues, genes);
  // copy share link
  await page.getByTestId(SHARE_BUTTON).click();
}

// eslint-disable-next-line sonarjs/cognitive-complexity
async function verifyShareLink(
  page: Page,
  version: string,
  tissues?: string[],
  genes?: string[],
  datasets?: Array<any>,
  sexes?: Array<any>,
  diseases?: Array<any>,
  ethnicities?: string[],
  _compare?: string
) {
  let encodedLink = `${TEST_URL}/gene-expression?`;
  let isFirstParam = false;
  let param = "";
  let data;
  // copy link to clipboard
  await page.getByTestId("share-button").click();
  const clipboardText: string = await page.evaluate(
    "navigator.clipboard.readText()"
  );

  // split parameters
  const urlParams = new URLSearchParams(
    // (thuang): We only want the query params part of the URL, so we split by "?"
    decodeURIComponent(clipboardText.split("?")[1])
  );

  // compare
  if (_compare !== undefined) {
    param = "compare";
    isFirstParam = true;
    await verifyParameter(page, urlParams, param, [_compare]);
    encodedLink += encodeLink(param, _compare, isFirstParam);
  }

  //datasets
  if (datasets !== undefined) {
    param = "datasets";
    data = await verifyParameter(page, urlParams, param, datasets || []);
    encodedLink += encodeLink(param, data.toString());
  }

  //diseases
  if (diseases !== undefined) {
    param = "diseases";
    data = await verifyParameter(page, urlParams, param, diseases || []);
    encodedLink += encodeLink(param, data.toString());
  }

  // verify ethnicities
  if (ethnicities !== undefined) {
    param = "ethnicities";
    data = await verifyParameter(page, urlParams, param, ethnicities || []);
    encodedLink += encodeLink(param, data?.toString());
  }

  // verify sexes
  if (sexes !== undefined) {
    param = "sexes";
    data = await verifyParameter(page, urlParams, param, sexes || []);
    encodedLink += encodeLink(param, data.toString());
  }

  // verify tissues
  if (tissues !== undefined) {
    param = "tissues";
    data = await verifyParameter(page, urlParams, param, tissues || []);
    if (!isFirstParam) {
      encodedLink += encodeLink(param, data?.toString(), true);
    } else {
      encodedLink += encodeLink(param, data?.toString());
    }
  }

  // verify genes
  if (genes !== undefined) {
    param = "genes";
    data = await verifyParameter(page, urlParams, param, genes || []);
    encodedLink += encodeLink(param, data?.toString());
  }

  // version
  param = "ver";
  data = await verifyParameter(page, urlParams, param, [version]);
  encodedLink += encodeLink(param, data?.toString());

  // verify encoded link
  expect(clipboardText).toBe(encodedLink);
}

function skipFirefox(browserName: string) {
  test.skip(browserName === "firefox", "No Clipboard read permission");
}
async function verifyParameter(
  page: Page,
  urlParams: URLSearchParams,
  param: string,
  data: Array<any>
): Promise<string[]> {
  if (data.length > 0) {
    let paramValues: string[] = [];
    // extact array od id attribute for {id: sting, text: string} data
    const ids = data.map((obj) => obj.id);
    switch (param) {
      case "datasets":
        paramValues = urlParams.get("datasets")?.split(",") || [];
        // verify datasets have been selected
        paramValues.forEach(async (_id: string) => {
          const item = data.find((item) => item.id === _id);
          if (item) {
            await expect(page.getByText(item.text)).toBeVisible();
          }
        });
        // verify query param values match
        expect(ids).toMatchObject(paramValues || {});
        break;
      case "tissues":
      case "genes":
      case "ethnicities":
      case "ver":
      case "compare":
        paramValues = urlParams.get(param)?.split(",") || [];
        //verify query param values match
        expect(data).toMatchObject(paramValues || {});
        break;
      default:
        paramValues = data.map((obj) => obj.id);
        expect(paramValues).toMatchObject(
          urlParams.get(param)?.split(",") || {}
        );
    }
    return paramValues;
  }
  return [];
}

function encodeLink(param: string, data: string, isFirstParam?: boolean) {
  let delimiter = "&";
  if (isFirstParam) {
    delimiter = "";
  }
  return `${delimiter}${param}=${encodeURIComponent(data)}`;
}
