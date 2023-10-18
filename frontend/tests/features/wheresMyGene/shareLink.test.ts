import { Page, expect } from "@playwright/test";
import { LATEST_SHARE_LINK_VERSION } from "src/views/WheresMyGeneV2/components/GeneSearchBar/components/ShareButton/utils";

import { TEST_URL } from "tests/common/constants";
import { test } from "tests/common/test";
import { tryUntil } from "tests/utils/helpers";

import { goToWMG, addTissuesAndGenes } from "tests/utils/wmgUtils";

const { describe } = test;

const SHARE_BUTTON = "share-button";

const BLOOD_TISSUE = {
  id: "UBERON:0000178",
  name: "blood",
};

const TISSUES = [BLOOD_TISSUE, { id: "UBERON:0002048", name: "lung" }];

const tissueNames = TISSUES.map((tissue) => tissue.name);
const tissueIds = TISSUES.map((tissue) => tissue.id);

const GENES = ["DPM1", "TNMD", "TSPAN6"];

const DATASETS = [
  {
    id: "d8da613f-e681-4c69-b463-e94f5e66847f",
    text: "A molecular single-cell lung atlas of lethal COVID-19",
  },
  {
    id: "de2c780c-1747-40bd-9ccf-9588ec186cee",
    text: "Immunophenotyping of COVID-19 and influenza highlights the role of type I interferons in development of severe COVID-19",
  },
];

const DISEASES = [{ id: "MONDO:0100096", text: "COVID-19" }];

const ETHNICITIES = ["unknown"];

const SEXES = [
  { id: "PATO:0000383", text: "female" },
  { id: "PATO:0000384", text: "male" },
];
const COMPARE = "disease";

const INCORRECT_CELL_TYPES = ["DOG", "CAT"];

const CELL_TYPES = ["CD4-positive helper T cell"];

const SHARE_LINK =
  `${TEST_URL}/gene-expression?` +
  // (thuang): Tissue params include a tissue id and a tissue name to test that we support both
  // (seve): cellType params includes incorrect cellTypes to test that we filter out invalid params
  `compare=disease&datasets=d8da613f-e681-4c69-b463-e94f5e66847f%2Cde2c780c-1747-40bd-9ccf-9588ec186cee&diseases=MONDO%3A0100096&ethnicities=unknown&sexes=PATO%3A0000383%2CPATO%3A0000384&tissues=blood%2CUBERON%3A0002048&genes=DPM1%2CTNMD%2CTSPAN&cellTypes=${INCORRECT_CELL_TYPES.concat(
    CELL_TYPES
  ).join("%2C")}&ver=2`;

describe("Share link tests", () => {
  test("Should share link with single tissue and single gene", async ({
    page,
    browserName,
  }) => {
    skipFirefox(browserName);
    const _genes = ["SCYL3"];

    //set up sate
    await setupStateAndCopyShareLink(page, [BLOOD_TISSUE.name], _genes);

    // verify link
    await verifyShareLink({
      page,
      linkVersion: LATEST_SHARE_LINK_VERSION,
      tissueIds: [BLOOD_TISSUE.id],
      genes: _genes,
    });
  });

  test("Should share link with multiple tissues and multiple genes", async ({
    page,
    browserName,
  }) => {
    skipFirefox(browserName);

    await setupStateAndCopyShareLink(page, tissueNames, GENES);

    // verify link
    await verifyShareLink({
      page,
      linkVersion: LATEST_SHARE_LINK_VERSION,
      tissueIds,
      genes: GENES,
    });
  });

  test("Should generate share link with correct format for all query param types", async ({
    page,
    browserName,
  }) => {
    skipFirefox(browserName);

    // prepare state
    await goToWMG(page, SHARE_LINK);

    // verify link parameters and app state
    await tryUntil(
      async () => {
        await verifyShareLink({
          page,
          linkVersion: LATEST_SHARE_LINK_VERSION,
          tissueIds,
          genes: GENES,
          datasets: DATASETS,
          sexes: SEXES,
          diseases: DISEASES,
          ethnicities: ETHNICITIES,
          compare: COMPARE,
          cellTypes: CELL_TYPES,
        });
      },
      { page }
    );
  });
});

async function setupStateAndCopyShareLink(
  page: Page,
  tissueNames: string[],
  genes: string[]
) {
  await goToWMG(page);
  await expect(page.getByTestId(SHARE_BUTTON)).toBeDisabled();

  // add tissues and genes
  await addTissuesAndGenes(page, tissueNames, genes);
  // copy share link
  await page.getByTestId(SHARE_BUTTON).click();
}

interface ExpectedParam {
  id: string;
  text: string;
}

async function verifyShareLink({
  page,
  linkVersion,
  tissueIds,
  genes,
  datasets,
  sexes,
  diseases,
  ethnicities,
  compare,
  cellTypes,
}: {
  page: Page;
  linkVersion: string;
  tissueIds?: string[];
  genes?: string[];
  datasets?: ExpectedParam[];
  sexes?: ExpectedParam[];
  diseases?: ExpectedParam[];
  ethnicities?: string[];
  compare?: string;
  cellTypes?: string[];
}) {
  let encodedLink = `${TEST_URL}/gene-expression?`;

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

  let isFirstParam = false;

  // compare
  if (compare !== undefined) {
    const param = "compare";
    isFirstParam = true;

    await verifyParameter(page, urlParams, param, [compare]);

    encodedLink += encodeLink(param, compare, isFirstParam);
  }

  // datasets
  if (datasets !== undefined) {
    const param = "datasets";

    const data = await verifyParameter(page, urlParams, param, datasets);

    encodedLink += encodeLink(param, String(data));
  }

  // diseases
  if (diseases !== undefined) {
    const param = "diseases";

    const data = await verifyParameter(page, urlParams, param, diseases);

    encodedLink += encodeLink(param, String(data));
  }

  // ethnicities
  if (ethnicities !== undefined) {
    const param = "ethnicities";

    const data = await verifyParameter(page, urlParams, param, ethnicities);

    encodedLink += encodeLink(param, String(data));
  }

  // sexes
  if (sexes !== undefined) {
    const param = "sexes";

    const data = await verifyParameter(page, urlParams, param, sexes);

    encodedLink += encodeLink(param, String(data));
  }

  // tissues
  if (tissueIds !== undefined) {
    const param = "tissues";

    const data = await verifyParameter(page, urlParams, param, tissueIds);

    if (!isFirstParam) {
      encodedLink += encodeLink(param, String(data), true);
    } else {
      encodedLink += encodeLink(param, String(data));
    }
  }

  // genes
  if (genes !== undefined) {
    const param = "genes";

    const data = await verifyParameter(page, urlParams, param, genes);

    encodedLink += encodeLink(param, String(data));
  }

  //cellTypes
  if (cellTypes !== undefined) {
    const param = "cellTypes";

    const data = await verifyParameter(page, urlParams, param, cellTypes);

    encodedLink += encodeLink(param, String(data));
  }

  // linkVersion
  const param = "ver";
  const data = await verifyParameter(page, urlParams, param, [linkVersion]);
  encodedLink += encodeLink(param, String(data));

  // encoded link
  expect(clipboardText).toBe(encodedLink);
}

function skipFirefox(browserName: string) {
  test.skip(browserName === "firefox", "No Clipboard read permission");
}
async function verifyParameter(
  page: Page,
  urlParams: URLSearchParams,
  param: string,
  expectedParams: Array<any>
): Promise<string[]> {
  if (!expectedParams.length) return [];

  // extract array of id attribute for {id: sting, text: string} expectedParams
  const expectedIds = expectedParams.map((expectedParam) => expectedParam.id);

  switch (param) {
    case "datasets": {
      const paramValues = getParamValues(param);

      // verify datasets have been selected
      paramValues.forEach(async (_id: string) => {
        const item = expectedParams.find(
          (expectedParam) => expectedParam.id === _id
        );

        if (item) {
          await expect(page.getByText(item.text)).toBeVisible();
        }
      });

      // verify query param values match
      expect(paramValues).toMatchObject(expectedIds);

      return paramValues;
    }
    case "tissues":
    case "genes":
    case "cellTypes":
    case "ethnicities":
    case "ver":
    case "compare": {
      const paramValues = getParamValues(param);

      // verify query param values match
      expect(paramValues).toMatchObject(expectedParams);

      return paramValues;
    }
    default: {
      const paramValues = getParamValues(param);

      expect(paramValues).toMatchObject(
        expectedParams.map((expectedParam) => expectedParam.id)
      );

      return paramValues;
    }
  }

  function getParamValues(param: string): string[] {
    return urlParams.get(param)?.split(",") || [];
  }
}

function encodeLink(param: string, data: string, isFirstParam?: boolean) {
  let delimiter = "&";
  if (isFirstParam) {
    delimiter = "";
  }
  return `${delimiter}${param}=${encodeURIComponent(data)}`;
}
