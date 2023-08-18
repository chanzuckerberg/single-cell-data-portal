import { test, Page, expect, Locator } from "@playwright/test";
import { ROUTES } from "src/common/constants/routes";
import {
  goToPage,
  takeSnapshotOfMetaTags,
  tryUntil,
} from "tests/utils/helpers";

import { TEST_URL } from "../../common/constants";

import {
  CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW,
  CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_HOVER_CONTAINER,
  CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_FULLSCREEN_BUTTON,
  CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_TOOLTIP,
} from "src/views/CellGuide/components/common/OntologyDagView/constants";

import { CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_CLICKABLE_TEXT_LABEL } from "src/views/CellGuide/components/common/OntologyDagView/components/Node/constants";

import { CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_RECT_OR_CIRCLE_PREFIX_ID } from "src/views/CellGuide/components/common/OntologyDagView/components/Node/components/RectOrCircle/constants";

import { CELL_GUIDE_CARD_SOURCE_DATA_TABLE } from "src/views/CellGuide/components/CellGuideCard/components/SourceDataTable/constants";

import { CELL_GUIDE_CARD_NAVIGATION_SIDEBAR } from "src/views/CellGuide/components/CellGuideCard/components/CellGuideCardSidebar/constants";

import {
  TISSUE_CARD_HEADER_NAME,
  TISSUE_CARD_HEADER_TAG,
  TISSUE_CARD_UBERON_DESCRIPTION,
} from "src/views/CellGuide/components/TissueCard/constants";

import {
  CELL_GUIDE_CARD_SEARCH_BAR,
  CELL_GUIDE_CARD_SEARCH_BAR_TEXT_INPUT,
} from "src/views/CellGuide/components/CellGuideCardSearchBar/constants";

import { LANDING_PAGE_HEADER } from "src/views/CellGuide/components/LandingPage/constants";

import {
  CELL_GUIDE_CARD_CL_DESCRIPTION,
  CELL_GUIDE_CARD_GPT_DESCRIPTION,
  CELL_GUIDE_CARD_GPT_TOOLTIP_LINK,
} from "src/views/CellGuide/components/CellGuideCard/components/Description/constants";

import {
  CELL_GUIDE_CARD_HEADER_NAME,
  CELL_GUIDE_CARD_HEADER_TAG,
  CELL_GUIDE_CARD_SYNONYMS,
} from "src/views/CellGuide/components/CellGuideCard/constants";
import {
  CELL_GUIDE_CARD_CANONICAL_MARKER_GENES_TABLE,
  CELL_GUIDE_CARD_CANONICAL_MARKER_GENES_TABLE_SELECTOR,
  CELL_GUIDE_CARD_ENRICHED_GENES_TABLE,
  CELL_GUIDE_CARD_ENRICHED_GENES_TABLE_SELECTOR,
  CELL_GUIDE_CARD_MARKER_GENES_TABLE_DROPDOWN_ORGAN,
  CELL_GUIDE_CARD_MARKER_GENES_TABLE_DROPDOWN_ORGANISM,
  EXPRESSION_SCORE_TOOLTIP_TEST_ID,
  MARKER_GENES_CANONICAL_TOOLTIP_TEST_ID,
  MARKER_GENES_COMPUTATIONAL_TOOLTIP_TEST_ID,
  PERCENT_OF_CELLS_TOOLTIP_TEST_ID,
} from "src/views/CellGuide/components/CellGuideCard/components/MarkerGeneTables/constants";

const { describe } = test;

const WAIT_FOR_TIMEOUT_MS = 30 * 1000;

const NEURON_CELL_TYPE_ID = "CL_0000540";
const T_CELL_CELL_TYPE_ID = "CL_0000084";
const LUNG_TISSUE_ID = "UBERON_0002048";

describe("Cell Guide", () => {
  describe("Landing Page", () => {
    test("All LandingPage components are present", async ({ page }) => {
      await goToPage(`${TEST_URL}${ROUTES.CELL_GUIDE}`, page);

      await takeSnapshotOfMetaTags("landing", page);

      await isElementVisible(page, LANDING_PAGE_HEADER);
      await isElementVisible(page, CELL_GUIDE_CARD_SEARCH_BAR);
    });

    test("Cell type search bar filters properly and links to a CellGuideCard", async ({
      page,
    }) => {
      await goToPage(`${TEST_URL}${ROUTES.CELL_GUIDE}`, page);

      const element = getSearchBarLocator(page);

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
      const firstOption = (await page.getByRole("option").all())[0];
      const firstOptionText = await firstOption?.textContent();
      expect(firstOptionText).toBe("neuron");

      await Promise.all([
        // check that the url has changed to the correct CellGuide card
        page.waitForURL(
          `${TEST_URL}${ROUTES.CELL_GUIDE}/${NEURON_CELL_TYPE_ID}`
        ),
        // click on first element in dropdown
        firstOption?.click(),
      ]);
    });
    test("Cell type search bar filters by synonyms properly and links to a CellGuideCard", async ({
      page,
    }) => {
      await goToPage(`${TEST_URL}${ROUTES.CELL_GUIDE}`, page);

      const element = getSearchBarLocator(page);

      await waitForElementAndClick(element);
      await waitForOptionsToLoad(page);
      // get number of elements with role option in dropdown
      const numOptionsBefore = await countLocator(page.getByRole("option"));
      // type in search bar
      await element.type("cardiomyocyte");
      // get number of elements with role option in dropdown
      const numOptionsAfter = await countLocator(page.getByRole("option"));
      // check that number of elements with role option in dropdown has decreased
      expect(numOptionsAfter).toBeLessThan(numOptionsBefore);
      // check that the first element in the dropdown is the one we searched for (cardiac muscle cell)
      const firstOption = (await page.getByRole("option").all())[0];
      const firstOptionText = await firstOption?.textContent();
      expect(firstOptionText).toBe("cardiac muscle cell");
      // click on first element in dropdown
      await firstOption?.click();
      // check that the url has changed to the correct CellGuide card
      await page.waitForURL(`${TEST_URL}${ROUTES.CELL_GUIDE}/CL_0000746`); // cardiac muscle cell
    });
    test("Cell type search bar keyboard input works properly", async ({
      page,
    }) => {
      await goToPage(`${TEST_URL}${ROUTES.CELL_GUIDE}`, page);

      const element = getSearchBarLocator(page);

      await waitForElementAndClick(element);

      await tryUntil(
        async () => {
          /**
           * (thuang): Clear the input field before typing in it.
           */
          await element.clear();
          await element.type("neuron");
          // input down arrow key
          await element.press("ArrowDown");
          const firstOption = page.getByRole("option").first();
          const firstOptionText = await firstOption?.textContent();
          expect(firstOptionText).toBe("neuron");

          // get css classes of first option
          const firstOptionClasses = await firstOption?.getAttribute("class");
          expect(firstOptionClasses).toContain("Mui-focused");
        },
        { page }
      );

      await Promise.all([
        // check that the url has changed to the correct CellGuide card after browser finishes navigating
        page.waitForURL(
          `${TEST_URL}${ROUTES.CELL_GUIDE}/${NEURON_CELL_TYPE_ID}`
        ), // input enter
        element.press("Enter"),
      ]);
    });
  });

  describe("CellGuide Card", () => {
    test("All CellGuide card components are present", async ({ page }) => {
      await goToPage(
        `${TEST_URL}${ROUTES.CELL_GUIDE}/${NEURON_CELL_TYPE_ID}`,
        page
      );

      await takeSnapshotOfMetaTags("cellType", page);

      await isElementVisible(page, CELL_GUIDE_CARD_HEADER_NAME);
      await isElementVisible(page, CELL_GUIDE_CARD_HEADER_TAG);
      await isElementVisible(page, CELL_GUIDE_CARD_CL_DESCRIPTION);
      await isElementVisible(page, CELL_GUIDE_CARD_GPT_DESCRIPTION);
      await isElementVisible(page, CELL_GUIDE_CARD_SYNONYMS);
      await isElementVisible(page, CELL_GUIDE_CARD_GPT_TOOLTIP_LINK);
      await isElementVisible(page, CELL_GUIDE_CARD_SEARCH_BAR);
      await isElementVisible(
        page,
        CELL_GUIDE_CARD_CANONICAL_MARKER_GENES_TABLE
      );
      await isElementVisible(page, CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW);
      const headerName = page.getByTestId(CELL_GUIDE_CARD_HEADER_NAME);
      const headerNameText = await headerName.textContent();
      expect(headerNameText).toBe("Neuron");
    });

    test("CellGuide card GPT description tooltip displays disclaimer", async ({
      page,
    }) => {
      await goToPage(
        `${TEST_URL}${ROUTES.CELL_GUIDE}/${NEURON_CELL_TYPE_ID}`,
        page
      );
      await isElementVisible(page, CELL_GUIDE_CARD_GPT_TOOLTIP_LINK);
      await page.getByTestId(CELL_GUIDE_CARD_GPT_TOOLTIP_LINK).hover();
      await checkTooltipContent(
        page,
        `This summary on "neuron" was generated with ChatGPT, powered by the GPT4 model. Keep in mind that ChatGPT may occasionally present information that is not entirely accurate. For transparency, the prompts used to generate this summary are shared below.`
      );
    });

    test("Cell type search bar filters properly and links to a CellGuideCard", async ({
      page,
    }) => {
      await goToPage(
        `${TEST_URL}${ROUTES.CELL_GUIDE}/${NEURON_CELL_TYPE_ID}`,
        page
      );
      const element = getSearchBarLocator(page);

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
      const firstOption = (await page.getByRole("option").all())[0];
      const firstOptionText = await firstOption?.textContent();
      expect(firstOptionText).toBe("acinar cell");

      await Promise.all([
        // check that the url has changed to the correct CellGuide card
        page.waitForURL(`${TEST_URL}${ROUTES.CELL_GUIDE}/CL_0000622`), // Acinar cell
        // click on first element in dropdown
        firstOption?.click(),
      ]);
    });

    describe("Canonical Marker Gene Table", () => {
      test("Canonical marker gene table is displayed with columns and at least one entry displayed", async ({
        page,
      }) => {
        await goToPage(
          `${TEST_URL}${ROUTES.CELL_GUIDE}/${NEURON_CELL_TYPE_ID}`,
          page
        );

        await tryUntil(
          async () => {
            // set canonical marker genes table as active
            await page
              .getByTestId(
                CELL_GUIDE_CARD_CANONICAL_MARKER_GENES_TABLE_SELECTOR
              )
              .click();

            const tableSelector = `[data-testid='${CELL_GUIDE_CARD_CANONICAL_MARKER_GENES_TABLE}']`;
            const columnHeaderElements = await page
              .locator(`${tableSelector} thead th`)
              .all();
            // get text content of each column header
            const columnHeaders = await Promise.all(
              columnHeaderElements.map(async (element) => {
                return await element.textContent();
              })
            );
            expect(columnHeaders).toEqual(["Symbol", "Name", "References"]);
            const rowElements = await page
              .locator(`${tableSelector} tbody tr`)
              .all();
            const rowCount = rowElements.length;
            expect(rowCount).toBeGreaterThan(1);
          },
          { page }
        );
      });

      test("Canonical marker gene table is updated by the tissue dropdown", async ({
        page,
      }) => {
        await goToPage(
          `${TEST_URL}${ROUTES.CELL_GUIDE}/${T_CELL_CELL_TYPE_ID}`,
          page
        );

        await tryUntil(
          async () => {
            // set canonical marker genes table as active
            await page
              .getByTestId(
                CELL_GUIDE_CARD_CANONICAL_MARKER_GENES_TABLE_SELECTOR
              )
              .click();

            const tableSelector = `[data-testid='${CELL_GUIDE_CARD_CANONICAL_MARKER_GENES_TABLE}']`;
            const rowElementsBefore = await page
              .locator(`${tableSelector} tbody tr`)
              .all();
            const rowCountBefore = rowElementsBefore.length;
            expect(rowCountBefore).toBeGreaterThan(1);

            const dropdown = page.getByTestId(
              CELL_GUIDE_CARD_MARKER_GENES_TABLE_DROPDOWN_ORGAN
            );
            await waitForElementAndClick(dropdown);
            await dropdown.press("ArrowDown");
            await dropdown.press("ArrowDown");
            await dropdown.press("ArrowDown"); // selects kidney
            await dropdown.press("Enter");

            const rowElementsAfter = await page
              .locator(`${tableSelector} tbody tr`)
              .all();
            const rowCountAfter = rowElementsAfter.length;
            expect(rowCountAfter).toBeGreaterThan(1);
            expect(rowCountAfter).not.toBe(rowCountBefore);
          },
          { page }
        );
      });
    });

    describe("Enriched Genes Table", () => {
      test("Enriched gene table is displayed with columns and at least one entry displayed", async ({
        page,
      }) => {
        await goToPage(
          `${TEST_URL}${ROUTES.CELL_GUIDE}/${T_CELL_CELL_TYPE_ID}`,
          page
        );

        const tableSelector = `[data-testid='${CELL_GUIDE_CARD_ENRICHED_GENES_TABLE}']`;

        await tryUntil(
          async () => {
            // set enriched marker genes table as active
            await page
              .getByTestId(CELL_GUIDE_CARD_ENRICHED_GENES_TABLE_SELECTOR)
              .click();

            await page
              .locator(tableSelector)
              .waitFor({ timeout: WAIT_FOR_TIMEOUT_MS });
          },
          { page }
        );

        const columnHeaderElements = await page
          .locator(`${tableSelector} thead th`)
          .all();

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
          .all();
        const rowCount = rowElements.length;
        expect(rowCount).toBeGreaterThan(1);
      });

      test("Enriched marker gene table is updated by the organism dropdown", async ({
        page,
      }) => {
        await goToPage(
          `${TEST_URL}${ROUTES.CELL_GUIDE}/${T_CELL_CELL_TYPE_ID}`,
          page
        );

        await tryUntil(
          async () => {
            // set enriched marker genes table as active
            await page
              .getByTestId(CELL_GUIDE_CARD_ENRICHED_GENES_TABLE_SELECTOR)
              .click();

            const tableSelector = `[data-testid='${CELL_GUIDE_CARD_ENRICHED_GENES_TABLE}']`;
            await page
              .locator(tableSelector)
              .waitFor({ timeout: WAIT_FOR_TIMEOUT_MS });

            const rowElementsBefore = await page
              .locator(`${tableSelector} tbody tr`)
              .all();
            const rowCountBefore = rowElementsBefore.length;
            expect(rowCountBefore).toBeGreaterThan(1);
            const firstRowContentBefore =
              await rowElementsBefore[0].textContent();

            const dropdown = page.getByTestId(
              CELL_GUIDE_CARD_MARKER_GENES_TABLE_DROPDOWN_ORGANISM
            );
            await waitForElementAndClick(dropdown);
            await dropdown.press("ArrowDown");
            await dropdown.press("Enter");

            const rowElementsAfter = await page
              .locator(`${tableSelector} tbody tr`)
              .all();
            const rowCountAfter = rowElementsAfter.length;
            expect(rowCountAfter).toBeGreaterThan(1);
            const firstRowContentAfter =
              await rowElementsAfter[0].textContent();
            expect(firstRowContentBefore).not.toBe(firstRowContentAfter);
          },
          { page }
        );
      });
      test("Enriched marker gene table is updated by the organ dropdown", async ({
        page,
      }) => {
        await goToPage(
          `${TEST_URL}${ROUTES.CELL_GUIDE}/${T_CELL_CELL_TYPE_ID}`,
          page
        );

        await tryUntil(
          async () => {
            // set enriched marker genes table as active
            await page
              .getByTestId(CELL_GUIDE_CARD_ENRICHED_GENES_TABLE_SELECTOR)
              .click();

            const tableSelector = `[data-testid='${CELL_GUIDE_CARD_ENRICHED_GENES_TABLE}']`;
            await page.locator(tableSelector).waitFor({ timeout: 5000 });

            const rowElementsBefore = await page
              .locator(`${tableSelector} tbody tr`)
              .all();
            const rowCountBefore = rowElementsBefore.length;
            expect(rowCountBefore).toBeGreaterThan(1);
            const firstRowContentBefore =
              await rowElementsBefore[0].textContent();

            const dropdown = page.getByTestId(
              CELL_GUIDE_CARD_MARKER_GENES_TABLE_DROPDOWN_ORGAN
            );
            await waitForElementAndClick(dropdown);
            await dropdown.press("ArrowDown");
            await dropdown.press("ArrowDown");
            await dropdown.press("Enter");

            const rowElementsAfter = await page
              .locator(`${tableSelector} tbody tr`)
              .all();
            const rowCountAfter = rowElementsAfter.length;
            expect(rowCountAfter).toBeGreaterThan(1);
            const firstRowContentAfter =
              await rowElementsAfter[0].textContent();
            expect(firstRowContentBefore).not.toBe(firstRowContentAfter);
          },
          { page }
        );
      });
    });

    describe("Source Data Table", () => {
      test("Source data table is displayed with columns and at least one entry displayed", async ({
        page,
      }) => {
        await goToPage(
          `${TEST_URL}${ROUTES.CELL_GUIDE}/${NEURON_CELL_TYPE_ID}`,
          page
        );

        const tableSelector = `[data-testid='${CELL_GUIDE_CARD_SOURCE_DATA_TABLE}']`;

        await tryUntil(
          async () => {
            const columnHeaderElements = await page
              .locator(`${tableSelector} thead th`)
              .all();

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
          },
          { page }
        );

        const rowElements = await page
          .locator(`${tableSelector} tbody tr`)
          .all();
        const rowCount = rowElements.length;
        expect(rowCount).toBeGreaterThan(1);
      });
    });

    describe("Ontology Viewer", () => {
      test("Clicking on a parent node expands and collapses its children", async ({
        page,
      }) => {
        await goToPage(
          `${TEST_URL}${ROUTES.CELL_GUIDE}/${NEURON_CELL_TYPE_ID}`,
          page
        );

        await page
          .getByTestId(CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW)
          .waitFor({ timeout: WAIT_FOR_TIMEOUT_MS });

        const nodesLocator = `[data-testid^='${CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_RECT_OR_CIRCLE_PREFIX_ID}']`;

        // collapse node's children
        const nodesBefore = await page.locator(nodesLocator).all();
        const numNodesBefore = nodesBefore.length;

        await tryUntil(
          async () => {
            const node = page.getByTestId(
              `${CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_RECT_OR_CIRCLE_PREFIX_ID}-CL:0000540__0-has-children-isTargetNode=true`
            );

            await waitForElementAndClick(node);

            const nodesAfter = await page.locator(nodesLocator).all();
            const numNodesAfter = nodesAfter.length;

            expect(numNodesBefore).toBeGreaterThan(numNodesAfter);

            // expand node's children
            await waitForElementAndClick(node);

            const nodesAfter2 = await page.locator(nodesLocator).all();
            const numNodesAfter2 = nodesAfter2.length;

            expect(numNodesAfter2).toBeGreaterThan(numNodesAfter);
          },
          { page }
        );
      });

      test("Clicking on a cell type label links to its CellGuide Card", async ({
        page,
      }) => {
        await goToPage(
          `${TEST_URL}${ROUTES.CELL_GUIDE}/${NEURON_CELL_TYPE_ID}`,
          page
        );
        await page
          .getByTestId(CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW)
          .waitFor({ timeout: WAIT_FOR_TIMEOUT_MS });
        const label = page.getByTestId(
          `${CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_CLICKABLE_TEXT_LABEL}-CL:0000878__4`
        );

        await Promise.all([
          page.waitForURL(`${TEST_URL}${ROUTES.CELL_GUIDE}/CL_0000878`),
          waitForElementAndClick(label),
        ]);

        // Check that the new node is highlighted green (isTargetNode=true)
        await page
          .getByTestId(
            `${CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_RECT_OR_CIRCLE_PREFIX_ID}-CL:0000878__4-has-children-isTargetNode=true`
          )
          .waitFor({ timeout: WAIT_FOR_TIMEOUT_MS });
      });

      test("Clicking on a collapsed node stub displays hidden cell types", async ({
        page,
      }) => {
        await goToPage(
          `${TEST_URL}${ROUTES.CELL_GUIDE}/${T_CELL_CELL_TYPE_ID}`,
          page
        );
        await page
          .getByTestId(CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW)
          .waitFor({ timeout: WAIT_FOR_TIMEOUT_MS });

        const nodesLocator = `[data-testid^='${CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_RECT_OR_CIRCLE_PREFIX_ID}']`;
        const nodesBefore = await page.locator(nodesLocator).all();
        const numNodesBefore = nodesBefore.length;

        // check that dummyChild is clickable
        const dummyChildLocator = page.getByTestId(
          `${CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_RECT_OR_CIRCLE_PREFIX_ID}-dummy-child-CL:0000842__0-has-children-isTargetNode=false`
        );
        await dummyChildLocator.click();

        const nodesAfter = await page.locator(nodesLocator).all();
        const numNodesAfter = nodesAfter.length;

        expect(numNodesAfter).toBeGreaterThan(numNodesBefore);
      });

      test("Clicking the full screen button maximizes the ontology viewer", async ({
        page,
      }) => {
        await goToPage(
          `${TEST_URL}${ROUTES.CELL_GUIDE}/${NEURON_CELL_TYPE_ID}`,
          page
        );
        const ontologyDagView = page.getByTestId(
          CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_HOVER_CONTAINER
        );
        await ontologyDagView.waitFor({ timeout: WAIT_FOR_TIMEOUT_MS });

        const ontologyDagViewSizeBefore = await ontologyDagView.boundingBox();

        await ontologyDagView.hover();

        const fullScreenButton = page.getByTestId(
          CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_FULLSCREEN_BUTTON
        );

        await tryUntil(
          async () => {
            await waitForElementAndClick(fullScreenButton);
          },
          { page }
        );

        const ontologyDagViewSizeAfter = await ontologyDagView.boundingBox();
        expect(ontologyDagViewSizeAfter?.width).toBeGreaterThan(
          ontologyDagViewSizeBefore?.width ?? 0
        );
      });

      test("Node tooltip displays on hover", async ({ page }) => {
        await goToPage(
          `${TEST_URL}${ROUTES.CELL_GUIDE}/${NEURON_CELL_TYPE_ID}`,
          page
        );
        await page
          .getByTestId(CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW)
          .waitFor({ timeout: WAIT_FOR_TIMEOUT_MS });

        const node = page.getByTestId(
          `${CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_RECT_OR_CIRCLE_PREFIX_ID}-CL:0000540__0-has-children-isTargetNode=true`
        );
        await node.hover();
        await isElementVisible(page, CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_TOOLTIP);
      });
    });
    describe("Tissue Card", () => {
      test("All tissue card components are present", async ({ page }) => {
        await goToPage(
          `${TEST_URL}${ROUTES.CELL_GUIDE}/tissues/${LUNG_TISSUE_ID}`,
          page
        );

        await takeSnapshotOfMetaTags("tissue", page);

        await isElementVisible(page, TISSUE_CARD_HEADER_NAME);
        await isElementVisible(page, TISSUE_CARD_HEADER_TAG);
        await isElementVisible(page, CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW);
        await isElementVisible(page, CELL_GUIDE_CARD_SEARCH_BAR);
        await isElementVisible(page, TISSUE_CARD_UBERON_DESCRIPTION);
        const headerName = page.getByTestId(TISSUE_CARD_HEADER_NAME);
        const headerNameText = await headerName.textContent();
        expect(headerNameText).toBe("Lung");
      });

      test("Clicking on a cell type label links to its CellGuide Card", async ({
        page,
      }) => {
        await goToPage(
          `${TEST_URL}${ROUTES.CELL_GUIDE}/tissues/${LUNG_TISSUE_ID}`,
          page
        );

        await page
          .getByTestId(CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW)
          .waitFor({ timeout: WAIT_FOR_TIMEOUT_MS });

        const label = page.getByTestId(
          `${CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_CLICKABLE_TEXT_LABEL}-CL:0000066__0`
        );

        await Promise.all([
          page.waitForURL(`${TEST_URL}${ROUTES.CELL_GUIDE}/CL_0000066`),
          waitForElementAndClick(label),
        ]);

        // Check that the new node is highlighted green (isTargetNode=true)
        await page
          .getByTestId(
            `${CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_RECT_OR_CIRCLE_PREFIX_ID}-CL:0000066__0-has-children-isTargetNode=true`
          )
          .waitFor({ timeout: WAIT_FOR_TIMEOUT_MS });
      });

      test("Node tooltip displays on hover", async ({ page }) => {
        await goToPage(
          `${TEST_URL}${ROUTES.CELL_GUIDE}/tissues/${LUNG_TISSUE_ID}`,
          page
        );
        await page
          .getByTestId(CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW)
          .waitFor({ timeout: WAIT_FOR_TIMEOUT_MS });

        const node = page.getByTestId(
          `${CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_RECT_OR_CIRCLE_PREFIX_ID}-CL:0000066__0-has-children-isTargetNode=false`
        );
        await node.hover();
        await isElementVisible(page, CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_TOOLTIP);
        const textContent = await page
          .getByTestId(CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_TOOLTIP)
          .textContent();
        expect(textContent).toContain("in lung");
      });
    });

    describe("CellGuideCard Sidebar", () => {
      test("Clicking on the navbar scrolls to the section", async ({
        page,
      }) => {
        await goToPage(
          `${TEST_URL}${ROUTES.CELL_GUIDE}/${NEURON_CELL_TYPE_ID}`,
          page
        );

        await tryUntil(
          async () => {
            const emptyState = await page
              .getByText(
                "marker genes for this cell type are unavailable at this time"
              )
              .all();
            expect(emptyState.length).toBe(0);
          },
          { page }
        );

        const navbar = page.getByTestId(CELL_GUIDE_CARD_NAVIGATION_SIDEBAR);

        // scroll to the bottom
        const section3 = page.getByTestId("section-3");
        await section3.scrollIntoViewIfNeeded();

        // check that source data is in viewport
        const sourceData = page.getByTestId(CELL_GUIDE_CARD_SOURCE_DATA_TABLE);
        await sourceData.waitFor({ timeout: WAIT_FOR_TIMEOUT_MS });
        expect(sourceData).toBeInViewport();

        // get the second navbar tab (ontology) and click to scroll
        const elements = await navbar
          .locator(".MuiButtonBase-root.MuiTab-root")
          .all();
        const tab = elements[1];
        await tab.click();

        // check that ontology is in viewport
        const ontologyView = page.getByTestId(
          CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW
        );
        await ontologyView.waitFor({ timeout: WAIT_FOR_TIMEOUT_MS });
        expect(ontologyView).toBeInViewport();

        // check that source data is not in viewport
        expect(sourceData).not.toBeInViewport();
      });
    });
  });
  describe("Marker Genes table shows tooltips", () => {
    const computationalTabButton =
      "cell-guide-card-enriched-genes-table-selector";

    test("Marker Genes help tooltip - Canonical", async ({ page }) => {
      await goToPage(
        `${TEST_URL}${ROUTES.CELL_GUIDE}/${NEURON_CELL_TYPE_ID}`,
        page
      );

      await isElementVisible(page, MARKER_GENES_CANONICAL_TOOLTIP_TEST_ID);
      await page.getByTestId(MARKER_GENES_CANONICAL_TOOLTIP_TEST_ID).hover();

      await checkTooltipContent(
        page,
        "Canonical marker genes and associated publications were derived from the Anatomical Structures, Cell Types and Biomarkers (ASCT+B) tables from the 5th Human Reference Atlas release (July 2023). The tables are authored and reviewed by an international team of anatomists, pathologists, physicians, and other experts."
      );
    });

    test("Marker Genes help tooltip - Computational", async ({ page }) => {
      await goToPage(
        `${TEST_URL}${ROUTES.CELL_GUIDE}/${NEURON_CELL_TYPE_ID}`,
        page
      );

      await page.getByTestId(computationalTabButton).click();

      await isElementVisible(page, MARKER_GENES_COMPUTATIONAL_TOOLTIP_TEST_ID);
      await page
        .getByTestId(MARKER_GENES_COMPUTATIONAL_TOOLTIP_TEST_ID)
        .hover();

      await checkTooltipContent(
        page,
        "Computational marker genes are derived from the CELLxGENE Census. They are computed utilizing the same methodology as featured in our Find Marker Genes feature from the Gene Expression application."
      );
    });

    test("Expression Score tooltip", async ({ page }) => {
      await goToPage(
        `${TEST_URL}${ROUTES.CELL_GUIDE}/${NEURON_CELL_TYPE_ID}`,
        page
      );

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
      await goToPage(
        `${TEST_URL}${ROUTES.CELL_GUIDE}/${NEURON_CELL_TYPE_ID}`,
        page
      );

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
  await tooltipLocator.waitFor({ timeout: WAIT_FOR_TIMEOUT_MS });
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
      await element.waitFor({ timeout: WAIT_FOR_TIMEOUT_MS });
      const isVisible = await element.isVisible();
      expect(isVisible).toBe(true);
    },
    { page }
  );
}

async function waitForElementAndClick(locator: Locator) {
  await locator.waitFor({ timeout: WAIT_FOR_TIMEOUT_MS });
  await locator.click();
}

async function countLocator(locator: Locator) {
  return (await locator.all()).length;
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

function getSearchBarLocator(page: Page) {
  return (
    page
      .getByTestId(CELL_GUIDE_CARD_SEARCH_BAR_TEXT_INPUT)
      /**
       * data-testid `CELL_GUIDE_CARD_SEARCH_BAR_TEXT_INPUT` is on a div instead of the <input />
       */
      .locator("input")
  );
}
