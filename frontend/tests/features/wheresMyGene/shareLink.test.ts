/* eslint-disable @typescript-eslint/no-explicit-any */
import { Page, expect, test } from "@playwright/test";
import { LATEST_SHARE_LINK_VERSION } from "src/views/WheresMyGene/components/GeneSearchBar/components/ShareButton/utils";
import { TEST_URL } from "tests/common/constants";
import { goToWMG, addTissuesAndGenes } from "tests/utils/wmgUtils";

const { describe } = test;

const tissues = ["blood", "lung", "liver"];
const genes = ["SCYL3", "TSPAN6", "TNMD"];
const dataSets = [
  {
    id: "c874f155-9bf9-4928-b821-f52c876b3e48",
    text: "49 years old male - Fresh PBMCs (1 day post-intubation)",
  },
  {
    id: "db59611b-42de-4035-93aa-1ed39f38b467",
    text: "49 years old male - Fresh PBMCs (2 day post-intubation)",
  },
  {
    id: "eeacb0c1-2217-4cf6-b8ce-1f0fedf1b569",
    text: "49 years old male - Fresh PBMCs (3 day post-intubation)",
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
const diseases = [{ id: "MONDO%3A0100096", text: "COVID-19" }];
const ethnicities = ["unknown"];
const sexes = [
  { id: "PATO%3A0000383", text: "female" },
  { id: "PATO%3A0000384", text: "male" },
];
const COMPARE = "disease";
const VERSION = "2";
const initialState =
  "https://localhost:3000/gene-expression?datasets=c874f155-9bf9-4928-b821-f52c876b3e48%2Cdb59611b-42de-4035-93aa-1ed39f38b467%2Ceeacb0c1-2217-4cf6-b8ce-1f0fedf1b569%2C881fe679-c6e0-45a3-9427-c4e81be6921f%2Cea786a06-5855-48b7-80d7-0313a21a2044%2C456e8b9b-f872-488b-871d-94534090a865&diseases=MONDO%3A0100096&ethnicities=unknown&sexes=PATO%3A0000383%2CPATO%3A0000384&tissues=blood%2Clung&genes=DPM1%2CTNMD%2CTSPAN6&ver=2&compare=disease";
describe("Share link tests", () => {
  test.only("Should share link with single tissue and single gene", async ({
    page,
    browserName,
  }) => {
    skipFirefox(browserName);

    await setupStateAndVerifyShareLink(page, ["blood"], ["SCYL3"]);
  });

  test("Should share link with multiple tissues and multiple genes", async ({
    page,
    browserName,
  }) => {
    skipFirefox(browserName);

    await setupStateAndVerifyShareLink(page, tissues, genes);
  });

  test.only("Should generate share link with correct format for all query param types", async ({
    page,
    browserName,
  }) => {
    skipFirefox(browserName);

    // prepare state
    await goToWMG(page, initialState);

    // verify link parameters and app state
    await verifyShareLink(
      page,
      tissues,
      genes,
      dataSets,
      sexes,
      diseases,
      ethnicities,
      VERSION,
      COMPARE
    );

    // copy share link
    const clipboardText: string = await page.evaluate(
      "navigator.clipboard.readText()"
    );

    await goToWMG(page, clipboardText);

    tissues.forEach(async (tissue) => {
      // selected tissue should be visible
      await expect(
        page.getByTestId(`cell-type-labels-${tissue}`)
      ).toBeVisible();
    });

    genes.forEach(async (gene) => {
      // selected gene should be visible
      expect(await page.getByTestId(`gene-name-${gene}`).textContent()).toBe(
        gene
      );
    });
  });
});

async function setupStateAndCopyShareLink(
  page: Page,
  tissues: string[],
  genes: string[]
) {
  await goToWMG(page);
  await expect(page.getByTestId("share-button")).toBeDisabled();

  // add tissues and genes
  await addTissuesAndGenes(page, tissues, genes);
  // copy share link
  await page.getByTestId("share-button").click();
}

async function setupStateAndVerifyShareLink(
  page: Page,
  tissues: string[],
  genes: string[]
) {
  await setupStateAndCopyShareLink(page, tissues, genes);

  // verify link
  await verifyShareLink(page, tissues, genes);
}

async function verifyShareLink(
  page: Page,
  tissues?: string[],
  genes?: string[],
  datasets?: Array<any>,
  sexes?: Array<any>,
  diseases?: Array<any>,
  ethnicities?: string[],
  ver?: string,
  _compare?: string
) {
  let encodedLink = `${TEST_URL}/gene-expression?`;

  // copy link to clipboard
  const clipboardText: string = await page.evaluate(
    "navigator.clipboard.readText()"
  );

  // split parameters
  const urlParams = new URLSearchParams(
    // (thuang): We only want the query params part of the URL, so we split by "?"
    decodeURIComponent(clipboardText.split("?")[1])
  );

  //datasets
  encodedLink += await verifyComplexParameter(
    page,
    urlParams,
    "datasets",
    datasets || []
  );
  //diseases
  encodedLink += await verifyComplexParameter(
    page,
    urlParams,
    "diseases",
    diseases || []
  );

  // verify ethnicities
  if (ethnicities !== undefined) {
    expect(ethnicities).toMatchObject(
      urlParams.get("ethnicities")?.split(",") || {}
    );
    encodedLink += `&ethnicities=${encodeURIComponent(ethnicities.toString())}`;
  }
  // verify sexes
  encodedLink += await verifyComplexParameter(
    page,
    urlParams,
    "sexes",
    sexes || []
  );
  // verify tissues
  if (tissues !== undefined) {
    expect(tissues).toMatchObject(urlParams.get("tissues")?.split(",") || {});
    encodedLink += `tissues=${encodeURIComponent(tissues.toString())}`;
  }

  // verify genes
  if (genes !== undefined) {
    expect(genes).toMatchObject(urlParams.get("genes")?.split(",") || {});
    encodedLink += `&genes=${encodeURIComponent(genes.toString())}`;
  }

  // compare
  if (_compare !== undefined) {
    encodedLink += verifySimpleParameter(urlParams, "compare", _compare);
  }


  // version
  const version = ver !== undefined ? ver : LATEST_SHARE_LINK_VERSION;
  encodedLink += verifySimpleParameter(urlParams, "ver", version);

  // verify encoded link
  expect(clipboardText).toBe(encodedLink);
}

function skipFirefox(browserName: string) {
  test.skip(browserName === "firefox", "No Clipboard read permission");
}
async function verifyComplexParameter(
  page: Page,
  urlParams: URLSearchParams,
  param: string,
  data: Array<any>
): Promise<string> {
  if (data.length > 0) {
    expect(Object.keys(data)).toMatchObject(
      urlParams.get(param)?.split(",") || {}
    );
    const ids: string[] = urlParams.get("datasets")?.split(",") || [];
    ids.forEach(async (_id: string) => {
      const item = data.find((item) => item.id === _id);
      if (item) {
        expect(page.getByText(item.text)).toBeVisible();
      }
    });
    const delimiter = "" ? param === "tissues" : "&";
    return `${delimiter}${param}=${encodeURIComponent(ids.toString())}`;
  }
  return "";
}
function verifySimpleParameter(
  urlParams: URLSearchParams,
  param: string,
  data: string
) {
  if (param !== undefined && param !== "" && data !== undefined) {
    console.log(param);
    console.log(urlParams.get(param)?.split(","));
    expect([data]).toMatchObject(
      urlParams.get(param)?.split(",") || {}
    );
    return `&${param}=${encodeURIComponent(data)}`;
  }
  return "";
}
