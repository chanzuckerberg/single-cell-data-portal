import { test, Page, expect, Locator } from "@playwright/test";
import { ROUTES } from "src/common/constants/routes";
import { goToPage, tryUntil } from "tests/utils/helpers";
import { TEST_URL } from "../../common/constants";
import { LANDING_PAGE_HEADER } from "src/views/CellCards/components/LandingPage";
import {
  CELL_CARD_SEARCH_BAR,
  CELL_CARD_SEARCH_BAR_TEXT_INPUT,
} from "src/views/CellCards/components/CellCardSearchBar";
import {
  CELL_CARD_CL_DESCRIPTION,
  CELL_CARD_GPT_DESCRIPTION,
  CELL_CARD_GPT_TOOLTIP_LINK,
} from "src/views/CellCards/components/CellCard/components/Description";
import {
  CELL_CARD_HEADER_NAME,
  CELL_CARD_HEADER_TAG,
} from "src/views/CellCards/components/CellCard";
import {
  CELL_CARD_CANONICAL_MARKER_GENES_TABLE,
  CELL_CARD_ENRICHED_GENES_TABLE,
  CELL_CARD_MARKER_GENES_TABLE_DROPDOWN_ORGAN,
  CELL_CARD_MARKER_GENES_TABLE_DROPDOWN_ORGANISM,
  CELL_CARD_CANONICAL_MARKER_GENES_TABLE_SELECTOR,
  CELL_CARD_ENRICHED_GENES_TABLE_SELECTOR,
  MARKER_GENES_TOOLTIP_TEST_ID,
} from "src/views/CellCards/components/CellCard/components/MarkerGeneTables";
import {
  CELL_CARD_ONTOLOGY_DAG_VIEW,
  CELL_CARD_ONTOLOGY_DAG_VIEW_HOVER_CONTAINER,
  CELL_CARD_ONTOLOGY_DAG_VIEW_FULLSCREEN_BUTTON,
  CELL_CARD_ONTOLOGY_DAG_VIEW_TOOLTIP,
} from "src/views/CellCards/components/common/OntologyDagView";
import { CELL_CARD_ONTOLOGY_DAG_VIEW_CLICKABLE_TEXT_LABEL } from "src/views/CellCards/components/common/OntologyDagView/components/Node";
import { CELL_CARD_ONTOLOGY_DAG_VIEW_RECT_OR_CIRCLE_PREFIX_ID } from "src/views/CellCards/components/common/OntologyDagView/components/Node/components/RectOrCircle";
import { CELL_CARD_SOURCE_DATA_TABLE } from "src/views/CellCards/components/CellCard/components/SourceDataTable";
import { CELL_CARD_NAVIGATION_SIDEBAR } from "src/views/CellCards/components/CellCard/components/CellCardSidebar";
import {
  TISSUE_CARD_HEADER_NAME,
  TISSUE_CARD_HEADER_TAG,
} from "src/views/CellCards/components/TissueCard";
import {
  EXPRESSION_SCORE_TOOLTIP_TEST_ID,
  PERCENT_OF_CELLS_TOOLTIP_TEST_ID,
} from "src/views/CellCards/components/CellCard/components/MarkerGeneTables";

const { describe } = test;

describe("Cell Cards", () => {
  describe("Landing Page", () => {
    test("All LandingPage components are present", async ({ page }) => {
      await goToPage(`${TEST_URL}${ROUTES.CELL_CARDS}`, page);
      await isElementVisible(page, LANDING_PAGE_HEADER);
      await isElementVisible(page, CELL_CARD_SEARCH_BAR);
    });
    test("Cell type search bar filters properly and links to a CellCard", async ({
      page,
    }) => {
      await goToPage(`${TEST_URL}${ROUTES.CELL_CARDS}`, page);
      const element = page.getByTestId(CELL_CARD_SEARCH_BAR_TEXT_INPUT);
      await waitForElementAndClick(element);
      await waitForOptionsToLoad(page);
      // get number of elements with role option in dropdown
      const numOptionsBefore = await countLocator(page.getByRole("option"));
      // type in search bar
      await element.type("neuron");
      // get number of elements with role option in dropdown
      const numOptionsAfter = await countLocator(page.getByRole("option"));
      // check that number of elements with role option in dropdown has decreased
      expect(numOptionsAfter).toBeLessThan(numOptionsBefore);
      // check that the first element in the dropdown is the one we searched for
      const firstOption = (await page.getByRole("option").elementHandles())[0];
      const firstOptionText = await firstOption?.textContent();
      expect(firstOptionText).toBe("neuron");
      // click on first element in dropdown
      await firstOption?.click();
      // check that the url has changed to the correct cell card
      await page.waitForURL(`${TEST_URL}${ROUTES.CELL_CARDS}/CL_0000540`); // Neuron
    });
    test("Cell type search bar keyboard input works properly", async ({
      page,
    }) => {
      await goToPage(`${TEST_URL}${ROUTES.CELL_CARDS}`, page);
      const element = page.getByTestId(CELL_CARD_SEARCH_BAR_TEXT_INPUT);
      await waitForElementAndClick(element);
      await element.type("neuron");
      // input down arrow key
      await element.press("ArrowDown");
      // input enter
      await element.press("Enter");
      // check that the url has changed to the correct cell card after browser finishes navigating
      await page.waitForURL(`${TEST_URL}${ROUTES.CELL_CARDS}/CL_0000540`); // Neuron
    });
  });

  describe("Cell Card", () => {
    test("All cell card components are present", async ({ page }) => {
      await goToPage(`${TEST_URL}${ROUTES.CELL_CARDS}/CL_0000540`, page); // Neuron
      await isElementVisible(page, CELL_CARD_HEADER_NAME);
      await isElementVisible(page, CELL_CARD_HEADER_TAG);
      await isElementVisible(page, CELL_CARD_CL_DESCRIPTION);
      await isElementVisible(page, CELL_CARD_GPT_DESCRIPTION);
      await isElementVisible(page, CELL_CARD_GPT_TOOLTIP_LINK);
      await isElementVisible(page, CELL_CARD_SEARCH_BAR);
      await isElementVisible(page, CELL_CARD_CANONICAL_MARKER_GENES_TABLE);
      await isElementVisible(page, CELL_CARD_ONTOLOGY_DAG_VIEW);
      const headerName = page.getByTestId(CELL_CARD_HEADER_NAME);
      const headerNameText = await headerName.textContent();
      expect(headerNameText).toBe("Neuron");
    });
    test("Cell card GPT description tooltip displays disclaimer", async ({
      page,
    }) => {
      await goToPage(`${TEST_URL}${ROUTES.CELL_CARDS}/CL_0000540`, page); // Neuron
      await isElementVisible(page, CELL_CARD_GPT_TOOLTIP_LINK);
      await page.getByTestId(CELL_CARD_GPT_TOOLTIP_LINK).hover();
      await checkTooltipContent(
        page,
        `This summary on "neuron" was generated with ChatGPT, powered by the GPT4 model. Keep in mind that ChatGPT may occasionally present information that is not entirely accurate. For transparency, the prompts used to generate this summary are shared below.`
      );
    });
    test("Cell type search bar filters properly and links to a CellCard", async ({
      page,
    }) => {
      await goToPage(`${TEST_URL}${ROUTES.CELL_CARDS}/CL_0000540`, page); // Neuron
      const element = page.getByTestId(CELL_CARD_SEARCH_BAR_TEXT_INPUT);
      await waitForElementAndClick(element);
      await waitForOptionsToLoad(page);
      // get number of elements with role option in dropdown
      const numOptionsBefore = await countLocator(page.getByRole("option"));
      // type in search bar
      await element.type("acinar cell");
      // get number of elements with role option in dropdown
      const numOptionsAfter = await countLocator(page.getByRole("option"));
      // check that number of elements with role option in dropdown has decreased
      expect(numOptionsAfter).toBeLessThan(numOptionsBefore);
      // check that the first element in the dropdown is the one we searched for
      const firstOption = (await page.getByRole("option").elementHandles())[0];
      const firstOptionText = await firstOption?.textContent();
      expect(firstOptionText).toBe("acinar cell");
      // click on first element in dropdown
      await firstOption?.click();
      // check that the url has changed to the correct cell card
      await page.waitForURL(`${TEST_URL}${ROUTES.CELL_CARDS}/CL_0000622`); // Acinar cell
    });
    describe("Canonical Marker Gene Table", () => {
      test("Canonical marker gene table is displayed with columns and at least one entry displayed", async ({
        page,
      }) => {
        await goToPage(`${TEST_URL}${ROUTES.CELL_CARDS}/CL_0000540`, page); // Neuron
        // set canonical marker genes table as active
        await page
          .getByTestId(CELL_CARD_CANONICAL_MARKER_GENES_TABLE_SELECTOR)
          .click();

        const tableSelector = `[data-testid='${CELL_CARD_CANONICAL_MARKER_GENES_TABLE}']`;
        const columnHeaderElements = await page
          .locator(`${tableSelector} thead th`)
          .elementHandles();
        // get text content of each column header
        const columnHeaders = await Promise.all(
          columnHeaderElements.map(async (element) => {
            return await element.textContent();
          })
        );
        expect(columnHeaders).toEqual(["Symbol", "Name", "References"]);
        const rowElements = await page
          .locator(`${tableSelector} tbody tr`)
          .elementHandles();
        const rowCount = rowElements.length;
        expect(rowCount).toBeGreaterThan(1);
      });
      test("Canonical marker gene table is updated by the tissue dropdown", async ({
        page,
      }) => {
        await goToPage(`${TEST_URL}${ROUTES.CELL_CARDS}/CL_0000084`, page); // T cell
        // set canonical marker genes table as active
        await page
          .getByTestId(CELL_CARD_CANONICAL_MARKER_GENES_TABLE_SELECTOR)
          .click();

        const tableSelector = `[data-testid='${CELL_CARD_CANONICAL_MARKER_GENES_TABLE}']`;
        const rowElementsBefore = await page
          .locator(`${tableSelector} tbody tr`)
          .elementHandles();
        const rowCountBefore = rowElementsBefore.length;
        expect(rowCountBefore).toBeGreaterThan(1);

        const dropdown = page.getByTestId(
          CELL_CARD_MARKER_GENES_TABLE_DROPDOWN_ORGAN
        );
        await waitForElementAndClick(dropdown);
        await dropdown.press("ArrowDown");
        await dropdown.press("ArrowDown");
        await dropdown.press("ArrowDown"); // selects kidney
        await dropdown.press("Enter");

        const rowElementsAfter = await page
          .locator(`${tableSelector} tbody tr`)
          .elementHandles();
        const rowCountAfter = rowElementsAfter.length;
        expect(rowCountAfter).toBeGreaterThan(1);
        expect(rowCountAfter).not.toBe(rowCountBefore);
      });
    });
    describe("Enriched Genes Table", () => {
      test("Enriched gene table is displayed with columns and at least one entry displayed", async ({
        page,
      }) => {
        await goToPage(`${TEST_URL}${ROUTES.CELL_CARDS}/CL_0000084`, page); // T cell
        // set enriched marker genes table as active
        await page.getByTestId(CELL_CARD_ENRICHED_GENES_TABLE_SELECTOR).click();

        const tableSelector = `[data-testid='${CELL_CARD_ENRICHED_GENES_TABLE}']`;
        await page.locator(tableSelector).waitFor({ timeout: 5000 });

        const columnHeaderElements = await page
          .locator(`${tableSelector} thead th`)
          .elementHandles();
        // get text content of each column header
        const columnHeaders = await Promise.all(
          columnHeaderElements.map(async (element) => {
            return await element.textContent();
          })
        );
        expect(columnHeaders).toEqual([
          "Symbol",
          "Name",
          "Marker Score",
          "Expression Score",
          "% of Cells",
        ]);
        const rowElements = await page
          .locator(`${tableSelector} tbody tr`)
          .elementHandles();
        const rowCount = rowElements.length;
        expect(rowCount).toBeGreaterThan(1);
      });
      test("Enriched marker gene table is updated by the organism dropdown", async ({
        page,
      }) => {
        await goToPage(`${TEST_URL}${ROUTES.CELL_CARDS}/CL_0000084`, page); // T cell
        // set enriched marker genes table as active
        await page.getByTestId(CELL_CARD_ENRICHED_GENES_TABLE_SELECTOR).click();

        const tableSelector = `[data-testid='${CELL_CARD_ENRICHED_GENES_TABLE}']`;
        await page.locator(tableSelector).waitFor({ timeout: 5000 });

        const rowElementsBefore = await page
          .locator(`${tableSelector} tbody tr`)
          .elementHandles();
        const rowCountBefore = rowElementsBefore.length;
        expect(rowCountBefore).toBeGreaterThan(1);
        const firstRowContentBefore = await rowElementsBefore[0].textContent();

        const dropdown = page.getByTestId(
          CELL_CARD_MARKER_GENES_TABLE_DROPDOWN_ORGANISM
        );
        await waitForElementAndClick(dropdown);
        await dropdown.press("ArrowDown");
        await dropdown.press("Enter");

        const rowElementsAfter = await page
          .locator(`${tableSelector} tbody tr`)
          .elementHandles();
        const rowCountAfter = rowElementsAfter.length;
        expect(rowCountAfter).toBeGreaterThan(1);
        const firstRowContentAfter = await rowElementsAfter[0].textContent();
        expect(firstRowContentBefore).not.toBe(firstRowContentAfter);
      });
    });

    describe("Source Data Table", () => {
      test("Source data table is displayed with columns and at least one entry displayed", async ({
        page,
      }) => {
        await goToPage(`${TEST_URL}${ROUTES.CELL_CARDS}/CL_0000540`, page); // Neuron
        const tableSelector = `[data-testid='${CELL_CARD_SOURCE_DATA_TABLE}']`;

        const columnHeaderElements = await page
          .locator(`${tableSelector} thead th`)
          .elementHandles();
        // get text content of each column header
        const columnHeaders = await Promise.all(
          columnHeaderElements.map(async (element) => {
            return await element.textContent();
          })
        );
        expect(columnHeaders).toEqual([
          "Collection",
          "Publication",
          "Tissue",
          "Disease",
          "Organism",
        ]);
        const rowElements = await page
          .locator(`${tableSelector} tbody tr`)
          .elementHandles();
        const rowCount = rowElements.length;
        expect(rowCount).toBeGreaterThan(1);
      });
    });

    describe("Ontology Viewer", () => {
      test("Clicking on a parent node expands and collapses its children", async ({
        page,
      }) => {
        await goToPage(`${TEST_URL}${ROUTES.CELL_CARDS}/CL_0000540`, page); // Neuron
        await page
          .getByTestId(CELL_CARD_ONTOLOGY_DAG_VIEW)
          .waitFor({ timeout: 5000 });
        const nodesLocator = `[data-testid^='${CELL_CARD_ONTOLOGY_DAG_VIEW_RECT_OR_CIRCLE_PREFIX_ID}']`;
        // collapse node's children
        const nodesBefore = await page.locator(nodesLocator).elementHandles();
        const numNodesBefore = nodesBefore.length;

        const node = page.getByTestId(
          `${CELL_CARD_ONTOLOGY_DAG_VIEW_RECT_OR_CIRCLE_PREFIX_ID}-CL:0000540__0-has-children-isTargetNode=true`
        );
        await waitForElementAndClick(node);

        const nodesAfter = await page.locator(nodesLocator).elementHandles();
        const numNodesAfter = nodesAfter.length;
        expect(numNodesBefore).toBeGreaterThan(numNodesAfter);

        // expand node's children
        await waitForElementAndClick(node);
        const nodesAfter2 = await page.locator(nodesLocator).elementHandles();
        const numNodesAfter2 = nodesAfter2.length;
        expect(numNodesAfter2).toBeGreaterThan(numNodesAfter);
      });
      test("Clicking on a cell type label links to its Cell Card", async ({
        page,
      }) => {
        await goToPage(`${TEST_URL}${ROUTES.CELL_CARDS}/CL_0000540`, page); // Neuron
        await page
          .getByTestId(CELL_CARD_ONTOLOGY_DAG_VIEW)
          .waitFor({ timeout: 5000 });
        const label = page.getByTestId(
          `${CELL_CARD_ONTOLOGY_DAG_VIEW_CLICKABLE_TEXT_LABEL}-CL:0000878__4`
        );
        await waitForElementAndClick(label);
        await page.waitForURL(`${TEST_URL}${ROUTES.CELL_CARDS}/CL_0000878`);
        // Check that the new node is highlighted green (isTargetNode=true)
        await page
          .getByTestId(
            `${CELL_CARD_ONTOLOGY_DAG_VIEW_RECT_OR_CIRCLE_PREFIX_ID}-CL:0000878__4-has-children-isTargetNode=true`
          )
          .waitFor({ timeout: 5000 });
      });
      test("Clicking on a collapsed node stub displays hidden cell types", async ({
        page,
      }) => {
        await goToPage(`${TEST_URL}${ROUTES.CELL_CARDS}/CL_0000084`, page); // T cell
        await page
          .getByTestId(CELL_CARD_ONTOLOGY_DAG_VIEW)
          .waitFor({ timeout: 5000 });

        const nodesLocator = `[data-testid^='${CELL_CARD_ONTOLOGY_DAG_VIEW_RECT_OR_CIRCLE_PREFIX_ID}']`;
        const nodesBefore = await page.locator(nodesLocator).elementHandles();
        const numNodesBefore = nodesBefore.length;

        const dummyChildLocator = `[data-testid='${CELL_CARD_ONTOLOGY_DAG_VIEW_RECT_OR_CIRCLE_PREFIX_ID}-dummy-child-CL:0000842__0-has-children-isTargetNode=false']`;
        const dummyChild = (
          await page.locator(dummyChildLocator).elementHandles()
        )[0];
        // check that dummyChild is clickable
        const isVisible = await dummyChild.isVisible();
        expect(isVisible).toBe(true);
        await dummyChild.click();

        const nodesAfter = await page.locator(nodesLocator).elementHandles();
        const numNodesAfter = nodesAfter.length;
        expect(numNodesAfter).toBeGreaterThan(numNodesBefore);
      });
      test("Clicking the full screen button maximizes the ontology viewer", async ({
        page,
      }) => {
        await goToPage(`${TEST_URL}${ROUTES.CELL_CARDS}/CL_0000540`, page); // Neuron
        const ontologyDagView = page.getByTestId(
          CELL_CARD_ONTOLOGY_DAG_VIEW_HOVER_CONTAINER
        );
        await ontologyDagView.waitFor({ timeout: 5000 });
        const ontologyDagViewSizeBefore = await ontologyDagView.boundingBox();
        await ontologyDagView.hover();
        const fullScreenButton = page.getByTestId(
          CELL_CARD_ONTOLOGY_DAG_VIEW_FULLSCREEN_BUTTON
        );
        await waitForElementAndClick(fullScreenButton);

        const ontologyDagViewSizeAfter = await ontologyDagView.boundingBox();
        expect(ontologyDagViewSizeAfter?.width).toBeGreaterThan(
          ontologyDagViewSizeBefore?.width ?? 0
        );
      });
      test("Node tooltip displays on hover", async ({ page }) => {
        await goToPage(`${TEST_URL}${ROUTES.CELL_CARDS}/CL_0000540`, page); // Neuron
        await page
          .getByTestId(CELL_CARD_ONTOLOGY_DAG_VIEW)
          .waitFor({ timeout: 5000 });

        const node = page.getByTestId(
          `${CELL_CARD_ONTOLOGY_DAG_VIEW_RECT_OR_CIRCLE_PREFIX_ID}-CL:0000540__0-has-children-isTargetNode=true`
        );
        await node.hover();
        await isElementVisible(page, CELL_CARD_ONTOLOGY_DAG_VIEW_TOOLTIP);
      });
    });
    describe("Tissue Card", () => {
      test("All tissue card components are present", async ({ page }) => {
        await goToPage(
          `${TEST_URL}${ROUTES.CELL_CARDS}/tissues/UBERON_0002048`,
          page
        ); // Lung
        await isElementVisible(page, TISSUE_CARD_HEADER_NAME);
        await isElementVisible(page, TISSUE_CARD_HEADER_TAG);
        await isElementVisible(page, CELL_CARD_ONTOLOGY_DAG_VIEW);
        await isElementVisible(page, CELL_CARD_SEARCH_BAR);
        const headerName = page.getByTestId(TISSUE_CARD_HEADER_NAME);
        const headerNameText = await headerName.textContent();
        expect(headerNameText).toBe("Lung");
      });
      test("Clicking on a cell type label links to its Cell Card", async ({
        page,
      }) => {
        await goToPage(
          `${TEST_URL}${ROUTES.CELL_CARDS}/tissues/UBERON_0002048`,
          page
        ); // Lung
        await page
          .getByTestId(CELL_CARD_ONTOLOGY_DAG_VIEW)
          .waitFor({ timeout: 5000 });
        const label = page.getByTestId(
          `${CELL_CARD_ONTOLOGY_DAG_VIEW_CLICKABLE_TEXT_LABEL}-CL:0000066__0`
        );
        await waitForElementAndClick(label);
        await page.waitForURL(`${TEST_URL}${ROUTES.CELL_CARDS}/CL_0000066`);
        // Check that the new node is highlighted green (isTargetNode=true)
        await page
          .getByTestId(
            `${CELL_CARD_ONTOLOGY_DAG_VIEW_RECT_OR_CIRCLE_PREFIX_ID}-CL:0000066__0-has-children-isTargetNode=true`
          )
          .waitFor({ timeout: 5000 });
      });
      test("Node tooltip displays on hover", async ({ page }) => {
        await goToPage(
          `${TEST_URL}${ROUTES.CELL_CARDS}/tissues/UBERON_0002048`,
          page
        ); // Lung
        await page
          .getByTestId(CELL_CARD_ONTOLOGY_DAG_VIEW)
          .waitFor({ timeout: 5000 });

        const node = page.getByTestId(
          `${CELL_CARD_ONTOLOGY_DAG_VIEW_RECT_OR_CIRCLE_PREFIX_ID}-CL:0000066__0-has-children-isTargetNode=false`
        );
        await node.hover();
        await isElementVisible(page, CELL_CARD_ONTOLOGY_DAG_VIEW_TOOLTIP);
        const textContent = await page
          .getByTestId(CELL_CARD_ONTOLOGY_DAG_VIEW_TOOLTIP)
          .textContent();
        expect(textContent).toContain("in lung");
      });
    });
    describe("CellCard Sidebar", () => {
      test("Clicking on the navbar scrolls to the section", async ({
        page,
      }) => {
        await goToPage(`${TEST_URL}${ROUTES.CELL_CARDS}/CL_0000540`, page); // Neuron
        const navbar = page.getByTestId(CELL_CARD_NAVIGATION_SIDEBAR);

        // scroll to the bottom
        const section3 = page.getByTestId("section-3");
        await section3.scrollIntoViewIfNeeded();

        // check that source data is in viewport
        const sourceData = page.getByTestId(CELL_CARD_SOURCE_DATA_TABLE);
        await sourceData.waitFor({ timeout: 5000 });
        expect(sourceData).toBeInViewport();

        // get the second navbar tab (ontology) and click to scroll
        const elements = await navbar
          .locator(".MuiButtonBase-root.MuiTab-root")
          .elementHandles();
        const tab = elements[1];
        await tab.click();

        // check that ontology is in viewport
        const ontologyView = page.getByTestId(CELL_CARD_ONTOLOGY_DAG_VIEW);
        await ontologyView.waitFor({ timeout: 5000 });
        expect(ontologyView).toBeInViewport();

        // check that source data is not in viewport
        expect(sourceData).not.toBeInViewport();
      });
    });
  });
  describe("Marker Genes table shows tooltips", () => {
    const computationalTabButton = "cell-card-enriched-genes-table-selector";

    test("Marker Genes help tooltip - Canonical", async ({ page }) => {
      await goToPage(`${TEST_URL}${ROUTES.CELL_CARDS}/CL_0000540`, page); // Neuron

      await isElementVisible(page, MARKER_GENES_TOOLTIP_TEST_ID);
      await page.getByTestId(MARKER_GENES_TOOLTIP_TEST_ID).hover();

      await checkTooltipContent(
        page,
        "The below marker genes and associated publications were derived from the Anatomical Structures, Cell Types and Biomarkers (ASCT+B) tables. The tables are authored and reviewed by an international team of anatomists, pathologists, physicians, and other experts."
      );
    });

    test("Marker Genes help tooltip - Computational", async ({ page }) => {
      await goToPage(`${TEST_URL}${ROUTES.CELL_CARDS}/CL_0000540`, page); // Neuron

      await page.getByTestId(computationalTabButton).click();

      await isElementVisible(page, MARKER_GENES_TOOLTIP_TEST_ID);
      await page.getByTestId(MARKER_GENES_TOOLTIP_TEST_ID).hover();

      await checkTooltipContent(
        page,
        "The marker genes listed below are computationally derived from the CELLxGENE corpus. They are computed utilizing the same methodology as featured in our Find Marker Genes feature from the Gene Expression application."
      );
    });

    test("Expression Score tooltip", async ({ page }) => {
      await goToPage(`${TEST_URL}${ROUTES.CELL_CARDS}/CL_0000540`, page); // Neuron

      await page.getByTestId(computationalTabButton).click();

      // Check expression score tooltip
      await isElementVisible(page, EXPRESSION_SCORE_TOOLTIP_TEST_ID);
      await page.getByTestId(EXPRESSION_SCORE_TOOLTIP_TEST_ID).hover();

      await checkTooltipContent(
        page,
        "The expression score is the average rankit-normalized gene expression among cells in the cell type that have non-zero values."
      );
    });

    test("Percent of Cells tooltip", async ({ page }) => {
      // MARKER_GENES_TOOLTIP_TEST_ID
      await goToPage(`${TEST_URL}${ROUTES.CELL_CARDS}/CL_0000540`, page); // Neuron

      await page.getByTestId(computationalTabButton).click();

      // Check Percent of Cells tooltip
      await isElementVisible(page, PERCENT_OF_CELLS_TOOLTIP_TEST_ID);
      await page.getByTestId(PERCENT_OF_CELLS_TOOLTIP_TEST_ID).hover();

      await checkTooltipContent(
        page,
        "Percentage of cells expressing a gene in the cell type. These numbers are calculated after cells with low coverage and low expression values have been filtered out."
      );
    });
  });
});

async function checkTooltipContent(page: Page, text: string) {
  // check role tooltip is visible
  const tooltipLocator = page.getByRole("tooltip");
  await tooltipLocator.waitFor({ timeout: 5000 });
  const tooltipLocatorVisible = await tooltipLocator.isVisible();
  expect(tooltipLocatorVisible).toBe(true);

  // check that tooltip contains text
  const tooltipText = await tooltipLocator.textContent();
  expect(tooltipText).toContain(text);
}

async function isElementVisible(page: Page, testId: string) {
  await tryUntil(
    async () => {
      const element = page.getByTestId(testId);
      await element.waitFor({ timeout: 5000 });
      const isVisible = await element.isVisible();
      expect(isVisible).toBe(true);
    },
    { page }
  );
}

async function waitForElementAndClick(locator: Locator) {
  await locator.waitFor({ timeout: 5000 });
  await locator.click();
}

async function countLocator(locator: Locator) {
  return (await locator.elementHandles()).length;
}

async function waitForOptionsToLoad(page: Page) {
  await tryUntil(
    async () => {
      const numOptions = await countLocator(page.getByRole("option"));
      expect(numOptions).toBeGreaterThan(0);
    },
    { page }
  );
}
