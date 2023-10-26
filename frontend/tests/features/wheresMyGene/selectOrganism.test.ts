import { Page, expect } from "playwright/test";
import { ROUTES } from "src/common/constants/routes";
import {
  ADD_TISSUE_ID,
  CELL_TYPE_FILTER_TEST_ID,
  COMPARE_DROPDOWN_ID,
  DATASET_FILTER_TEST_ID,
  DISEASE_FILTER_TEST_ID,
  PUBLICATION_FILTER_TEST_ID,
  SELF_REPORTED_ETHNICITY_FILTER_TEST_ID,
  SEX_FILTER_TEST_ID,
  TEST_URL,
} from "tests/common/constants";
import { test } from "tests/common/test";
import { tryUntil } from "tests/utils/helpers";
import { WMG_WITH_SEEDED_GENES, goToWMG } from "tests/utils/wmgUtils";

const { describe } = test;

const TISSUE_NAME_TEST_ID = "tissue-name";
const MOUSE_NAME = "Mus musculus";

describe("Select organism", () => {
  test("Switch organism to `Mus musculus` with no customization", async ({
    page,
  }) => {
    const { tissueNameCount: beforeTissueNameCount } = await setup({ page });

    const targetOption = page.getByText(MOUSE_NAME);

    await tryUntil(
      async () => {
        await page.getByTestId("add-organism").click();
        await expect(targetOption).toBeVisible();
      },
      { page }
    );

    await targetOption.click();

    await tryUntil(
      async () => {
        expect(await page.getByTestId(TISSUE_NAME_TEST_ID).count()).not.toBe(
          beforeTissueNameCount
        );
      },
      { page }
    );
  });

  describe("Switch organism to `Mus musculus` with customizations", () => {
    test("With genes selected, it prompts confirm modal on change", async ({
      page,
    }) => {
      await verifyConfirmModal({ page, url: WMG_WITH_SEEDED_GENES.URL });
    });

    test("With dataset `22 integrated samples` selected, it prompts confirm modal on change", async ({
      page,
    }) => {
      await verifyConfirmModal({
        page,
        url:
          `${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}?` +
          "datasets=9d8e5dca-03a3-457d-b7fb-844c75735c83&ver=2",
        testId: DATASET_FILTER_TEST_ID,
        text: "22 integrated samples",
      });
    });

    test("With disease `normal` selected, it prompts confirm modal on change", async ({
      page,
    }) => {
      await verifyConfirmModal({
        page,
        url:
          `${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}?` +
          "diseases=PATO%3A0000461&ver=2",
        testId: DISEASE_FILTER_TEST_ID,
        text: "normal",
      });
    });

    test("With ethnicity `unknown` selected, it prompts confirm modal on change", async ({
      page,
    }) => {
      await verifyConfirmModal({
        page,
        url:
          `${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}?` +
          "ethnicities=unknown&ver=2",
        testId: SELF_REPORTED_ETHNICITY_FILTER_TEST_ID,
        text: "unknown",
      });
    });

    test("With publication `Ahern et al.` selected, it prompts confirm modal on change", async ({
      page,
    }) => {
      await verifyConfirmModal({
        page,
        url:
          `${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}?` +
          "publications=Ahern+et+al.+Cell+2022&ver=2",
        testId: PUBLICATION_FILTER_TEST_ID,
        text: "Ahern et al.",
      });
    });

    test("With sex `unknown` selected, it prompts confirm modal on change", async ({
      page,
    }) => {
      await verifyConfirmModal({
        page,
        url: `${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}?` + "sexes=unknown&ver=2",
        testId: SEX_FILTER_TEST_ID,
        text: "unknown",
      });
    });

    test("With tissue `lung` selected, it prompts confirm modal on change", async ({
      page,
    }) => {
      await verifyConfirmModal({
        page,
        url:
          `${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}?` +
          "tissues=UBERON%3A0002048&ver=2",
        testId: ADD_TISSUE_ID,
        text: "lung",
      });
    });

    test("With groupBy `disease` selected, it prompts confirm modal on change", async ({
      page,
    }) => {
      await verifyConfirmModal({
        page,
        url:
          `${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}?` +
          /**
           * (thuang): Somehow compare share link needs a gene to work
           */
          "compare=disease&genes=TSPAN6&ver=2",
        testId: COMPARE_DROPDOWN_ID,
        text: "Disease",
        /**
         * (thuang): Explicitly remove the seeded gene, so we can test just the groupBy
         */
        removeGene: "TSPAN6",
      });
    });

    test("With cellType `B cell` selected, it prompts confirm modal on change", async ({
      page,
    }) => {
      await verifyConfirmModal({
        page,
        url:
          `${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}?` +
          /**
           * (thuang): Somehow compare share link needs a gene to work
           */
          "cellTypes=B+cell&genes=TSPAN6&ver=2",
        testId: CELL_TYPE_FILTER_TEST_ID,
        text: "B cell",
        /**
         * (thuang): Explicitly remove the seeded gene, so we can test just the cellType
         */
        removeGene: "TSPAN6",
      });
    });
  });
});

async function verifyConfirmModal({
  page,
  url,
  testId,
  text,
  removeGene,
}: {
  page: Page;
  url?: string;
  testId?: string;
  text?: string;
  removeGene?: string;
}) {
  const { tissueNameCount: beforeTissueNameCount } = await setup({
    page,
    url,
  });

  if (testId && text) {
    await tryUntil(
      async () => {
        await expect(page.getByTestId(testId).getByText(text)).toBeVisible();
      },
      { page }
    );
  }

  if (removeGene) {
    await page.getByTestId("gene-name-" + removeGene).hover();
    await page.getByTestId("gene-delete-icon-" + removeGene).click();
  }

  const targetOption = page.getByText(MOUSE_NAME);

  await tryUntil(
    async () => {
      await page.getByTestId("add-organism").click();
      await expect(targetOption).toBeVisible();
    },
    { page }
  );

  await targetOption.click();

  await expect(page.getByRole("dialog")).toBeVisible();

  await page.getByText("Confirm").click();

  await tryUntil(
    async () => {
      expect(await page.getByTestId("tissue-name").count()).not.toBe(
        beforeTissueNameCount
      );
    },
    { page }
  );
}

async function setup({ page, url }: { page: Page; url?: string }): Promise<{
  tissueNameCount: number;
}> {
  await goToWMG(page, url);

  const tissueNameSelector = page.getByTestId(TISSUE_NAME_TEST_ID);

  await tryUntil(
    async () => {
      expect(await tissueNameSelector.count()).toBeGreaterThan(0);
    },
    { page }
  );

  return {
    tissueNameCount: await tissueNameSelector.count(),
  };
}
