import { Page, expect, Locator } from "@playwright/test";
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
  CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_DEACTIVATE_MARKER_GENE_MODE,
} from "src/views/CellGuide/components/common/OntologyDagView/constants";
import { CELL_GUIDE_ONTOLOGY_VIEW_LEGEND_TEST_ID } from "src/views/CellGuide/components/common/OntologyDagView/components/Legend/constants";
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
  CELL_GUIDE_CARD_SYNONYMS,
  CELL_GUIDE_CARD_VALIDATED_DESCRIPTION,
} from "src/views/CellGuide/components/CellGuideCard/components/Description/constants";

import {
  CELL_GUIDE_CARD_GLOBAL_ORGANISM_FILTER_DROPDOWN,
  CELL_GUIDE_CARD_GLOBAL_TISSUE_FILTER_DROPDOWN,
  CELL_GUIDE_CARD_HEADER_NAME,
  CELL_GUIDE_CARD_HEADER_TAG,
} from "src/views/CellGuide/components/CellGuideCard/constants";

import {
  CELL_GUIDE_CARD_CANONICAL_MARKER_GENES_TABLE,
  CELL_GUIDE_CARD_CANONICAL_MARKER_GENES_TABLE_SELECTOR,
  CELL_GUIDE_CARD_ENRICHED_GENES_TABLE,
  CELL_GUIDE_CARD_ENRICHED_GENES_TABLE_SELECTOR,
  EXPRESSION_SCORE_TOOLTIP_TEST_ID,
  MARKER_GENES_CANONICAL_TOOLTIP_TEST_ID,
  MARKER_GENES_COMPUTATIONAL_TOOLTIP_TEST_ID,
  PERCENT_OF_CELLS_TOOLTIP_TEST_ID,
  MARKER_GENES_TREE_ICON_BUTTON_TEST_ID,
} from "src/views/CellGuide/components/CellGuideCard/components/MarkerGeneTables/constants";
import {
  CELLGUIDE_INFO_SIDEBAR_TEST_ID,
  CELLGUIDE_VIEW_PAGE_SIDEBAR_BUTTON_TEST_ID,
} from "src/views/CellGuide/components/CellGuideInfoSideBar/constants";
import { test } from "tests/common/test";

const { describe } = test;

const WAIT_FOR_TIMEOUT_MS = 30 * 1000;

const NEURON_CELL_TYPE_ID = "CL_0000540";
const GLIOBLAST_CELL_TYPE_ID = "CL_0000030";
const T_CELL_CELL_TYPE_ID = "CL_0000084";
const LUNG_TISSUE_ID = "UBERON_0002048";
const ABNORMAL_CELL_TYPE_ID = "CL_0001061";

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
      await isElementVisible(page, CELL_GUIDE_CARD_ENRICHED_GENES_TABLE);
      await isElementVisible(page, CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW);
      const headerName = page.getByTestId(CELL_GUIDE_CARD_HEADER_NAME);
      const headerNameText = await headerName.textContent();
      expect(headerNameText).toBe("Neuron");
    });
    test("Glioblast CellGuide card is validated", async ({ page }) => {
      await goToPage(
        `${TEST_URL}${ROUTES.CELL_GUIDE}/${GLIOBLAST_CELL_TYPE_ID}`,
        page
      );
      await isElementVisible(page, CELL_GUIDE_CARD_VALIDATED_DESCRIPTION);
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
      /**
       * (thuang): Sometimes the dropdown options don't load, so add a tryUntil
       * to refresh the page and try again.
       */
      await tryUntil(
        async () => {
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
        },
        { page }
      );
    });

    describe("Canonical Marker Gene Table", () => {
      test("Canonical marker gene table is displayed with columns and at least one entry displayed", async ({
        page,
      }) => {
        await goToPage(
          `${TEST_URL}${ROUTES.CELL_GUIDE}/${NEURON_CELL_TYPE_ID}`,
          page
        );

        const tableSelector = `[data-testid='${CELL_GUIDE_CARD_CANONICAL_MARKER_GENES_TABLE}']`;
        await selectCanonicalMarkerGeneTable(page);

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
        const tableSelector = `[data-testid='${CELL_GUIDE_CARD_CANONICAL_MARKER_GENES_TABLE}']`;
        await selectCanonicalMarkerGeneTable(page);
        await tryUntil(
          async () => {
            const rowElementsBefore = await page
              .locator(`${tableSelector} tbody tr`)
              .all();
            const rowCountBefore = rowElementsBefore.length;
            expect(rowCountBefore).toBeGreaterThan(1);

            const dropdown = page.getByTestId(
              CELL_GUIDE_CARD_GLOBAL_TISSUE_FILTER_DROPDOWN
            );

            await tryUntil(
              async () => {
                await waitForElementAndClick(dropdown);
                await page.getByRole("option").getByText("brain").click();
                const rowElementsAfter = await page
                  .locator(`${tableSelector} tbody tr`)
                  .all();

                const rowCountAfter = rowElementsAfter.length;
                expect(rowCountAfter).toBeGreaterThan(1);
                expect(rowCountAfter).not.toBe(rowCountBefore);
              },
              { page }
            );
          },
          { page }
        );
      });

      test("Canonical marker gene table is updated by the organism dropdown", async ({
        page,
      }) => {
        await goToPage(
          `${TEST_URL}${ROUTES.CELL_GUIDE}/${T_CELL_CELL_TYPE_ID}`,
          page
        );
        await selectCanonicalMarkerGeneTable(page);
        const tableSelector = `[data-testid='${CELL_GUIDE_CARD_CANONICAL_MARKER_GENES_TABLE}']`;

        await tryUntil(
          async () => {
            const rowElementsBefore = await page
              .locator(`${tableSelector} tbody tr`)
              .all();
            const rowCountBefore = rowElementsBefore.length;
            expect(rowCountBefore).toBeGreaterThanOrEqual(1);

            const dropdown = page.getByTestId(
              CELL_GUIDE_CARD_GLOBAL_ORGANISM_FILTER_DROPDOWN
            );
            await waitForElementAndClick(dropdown);
            await page.getByRole("option").getByText("Mus musculus").click();

            const rowElementsAfter = await page
              .locator(`${tableSelector} tbody tr`)
              .all();
            const rowCountAfter = rowElementsAfter.length;
            expect(rowCountAfter).toBe(0);
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

        await selectComputationalMarkerGeneTable(page);

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

        await selectComputationalMarkerGeneTable(page);
        const tableSelector = `[data-testid='${CELL_GUIDE_CARD_ENRICHED_GENES_TABLE}']`;

        await tryUntil(
          async () => {
            const rowElementsBefore = await page
              .locator(`${tableSelector} tbody tr`)
              .all();
            const rowCountBefore = rowElementsBefore.length;
            expect(rowCountBefore).toBeGreaterThan(1);
            const firstRowContentBefore =
              await rowElementsBefore[0].textContent();

            const dropdown = page.getByTestId(
              CELL_GUIDE_CARD_GLOBAL_ORGANISM_FILTER_DROPDOWN
            );
            await waitForElementAndClick(dropdown);
            await page.getByRole("option").getByText("Mus musculus").click();

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

        await selectComputationalMarkerGeneTable(page);
        const tableSelector = `[data-testid='${CELL_GUIDE_CARD_ENRICHED_GENES_TABLE}']`;

        await tryUntil(
          async () => {
            const rowElementsBefore = await page
              .locator(`${tableSelector} tbody tr`)
              .all();
            const rowCountBefore = rowElementsBefore.length;
            expect(rowCountBefore).toBeGreaterThan(1);
            const firstRowContentBefore =
              await rowElementsBefore[0].textContent();

            const dropdown = page.getByTestId(
              CELL_GUIDE_CARD_GLOBAL_TISSUE_FILTER_DROPDOWN
            );
            await waitForElementAndClick(dropdown);
            await page.keyboard.type("colon");
            await page.keyboard.press("ArrowDown");
            await page.keyboard.press("Enter");

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
        /**
         * (thuang): Sometimes the table doesn't load, so add a tryUntil to refresh
         * the page and try again.
         */
        await tryUntil(
          async () => {
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
          },
          { page }
        );
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

        // Collapse node's children
        await tryUntil(
          async () => {
            const nodesBefore = await page.locator(nodesLocator).all();
            const numNodesBefore = nodesBefore.length;

            /**
             * (thuang): This is needed to ensure that we don't query the tree
             * before it's rendered
             */
            expect(numNodesBefore).toBeGreaterThan(0);

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

      test("Clicking on a cell type label opens its CellGuide info sidebar", async ({
        page,
      }) => {
        await goToPage(
          `${TEST_URL}${ROUTES.CELL_GUIDE}/${NEURON_CELL_TYPE_ID}`,
          page
        );
        await page
          .getByTestId(CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW)
          .waitFor({ timeout: WAIT_FOR_TIMEOUT_MS });
        await page
          .getByTestId(
            `${CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_CLICKABLE_TEXT_LABEL}-CL:0000878__4`
          )
          .click();

        await isElementVisible(page, CELLGUIDE_INFO_SIDEBAR_TEST_ID);
        await page
          .getByTestId(CELLGUIDE_VIEW_PAGE_SIDEBAR_BUTTON_TEST_ID)
          .click();

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
        await page.getByTestId("section-1").scrollIntoViewIfNeeded();
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

      test("Tree view exists for a cell type that is not a descendant of animal cell", async ({
        page,
      }) => {
        await goToPage(
          `${TEST_URL}${ROUTES.CELL_GUIDE}/${ABNORMAL_CELL_TYPE_ID}`,
          page
        );
        await page
          .getByTestId(CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW)
          .waitFor({ timeout: WAIT_FOR_TIMEOUT_MS });

        const node = page.getByTestId(
          `${CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_RECT_OR_CIRCLE_PREFIX_ID}-${ABNORMAL_CELL_TYPE_ID.replace(
            "_",
            ":"
          )}__0-has-children-isTargetNode=true`
        );
        await node.hover();
        await isElementVisible(page, CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_TOOLTIP);
        // check that the node `cell` (CL:0000000) is displayed
        await isElementVisible(
          page,
          `${CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_RECT_OR_CIRCLE_PREFIX_ID}-CL:0000000__0-has-children-isTargetNode=false`
        );
      });

      test("Clicking on a computational marker gene tree icon enters marker gene mode in a CellGuide Card", async ({
        page,
      }) => {
        await goToPage(
          `${TEST_URL}${ROUTES.CELL_GUIDE}/${NEURON_CELL_TYPE_ID}`,
          page
        );
        await page
          .getByTestId(CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW)
          .waitFor({ timeout: WAIT_FOR_TIMEOUT_MS });

        const tableSelector = `[data-testid='${CELL_GUIDE_CARD_ENRICHED_GENES_TABLE}']`;

        await selectComputationalMarkerGeneTable(page);

        const rowElements = await page
          .locator(`${tableSelector} tbody tr`)
          .all();
        const rowText = await rowElements[0].textContent();
        const geneSymbol = rowText?.split(" ").at(0);
        expect(geneSymbol).toBeDefined();
        const treeIcon = page.getByTestId(
          MARKER_GENES_TREE_ICON_BUTTON_TEST_ID(geneSymbol as string)
        );

        await rowElements[0].locator("td").nth(0).hover();
        await treeIcon.click();

        // check that the eyeClosed button is visible
        await isElementVisible(
          page,
          CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_DEACTIVATE_MARKER_GENE_MODE
        );

        // hover over the node
        const node = page.getByTestId(
          `${CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_RECT_OR_CIRCLE_PREFIX_ID}-CL:0000123__0-has-children-isTargetNode=false`
        );
        await node.hover();
        await isElementVisible(page, CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_TOOLTIP);

        // assert that the tooltip text contains the marker gene information
        const tooltipText = await page
          .getByTestId(CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_TOOLTIP)
          .textContent();

        expect(tooltipText).toContain(`${geneSymbol} stats`);

        const legendText = await page
          .getByTestId(CELL_GUIDE_ONTOLOGY_VIEW_LEGEND_TEST_ID)
          .textContent();
        expect(legendText).toContain("Marker Score");
        expect(legendText).toContain("Expressed in Cells(%)");

        // deactivate marker gene mode and check that the legend and tooltips reverted
        await page
          .getByTestId(
            CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_DEACTIVATE_MARKER_GENE_MODE
          )
          .click();

        const node2 = page.getByTestId(
          `${CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_RECT_OR_CIRCLE_PREFIX_ID}-CL:0000047__0-has-children-isTargetNode=false`
        );
        await node2.hover();
        await isElementVisible(page, CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_TOOLTIP);

        const newTooltipText = await page
          .getByTestId(CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_TOOLTIP)
          .textContent();

        expect(newTooltipText).not.toContain(`${geneSymbol} stats`);
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
          `${CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_CLICKABLE_TEXT_LABEL}-CL:1000271__0`
        );

        await Promise.all([
          page.waitForURL(`${TEST_URL}${ROUTES.CELL_GUIDE}/CL_1000271`),
          waitForElementAndClick(label),
        ]);

        // Check that the new node is highlighted green (isTargetNode=true)
        await page
          .getByTestId(
            `${CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_RECT_OR_CIRCLE_PREFIX_ID}-CL:1000271__0-no-children-isTargetNode=true`
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
          `${CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW_RECT_OR_CIRCLE_PREFIX_ID}-CL:1000271__0-no-children-isTargetNode=false`
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
        const sourceData = page.getByTestId(CELL_GUIDE_CARD_SOURCE_DATA_TABLE);
        const ontologyView = page.getByTestId(
          CELL_GUIDE_CARD_ONTOLOGY_DAG_VIEW
        );
        // scroll to the bottom

        const section3 = page.getByTestId("section-3");
        await section3.scrollIntoViewIfNeeded();
        await tryUntil(
          async () => {
            // check that source data is in viewport
            expect(sourceData).toBeInViewport();
            // 1 second in between retries or we hit the retry limit too fast
            await page.waitForTimeout(1000);
          },
          { page }
        );
        // get the second navbar tab (ontology) and click to scroll
        const elements = await navbar
          .locator(".MuiButtonBase-root.MuiTab-root")
          .all();
        const tab = elements[1];
        await tab.click();
        await tryUntil(
          async () => {
            // check that ontology is in viewport
            expect(ontologyView).toBeInViewport();
            // check that source data is not in viewport
            expect(sourceData).not.toBeInViewport();
            // 1 second in between retries or we hit the retry limit too fast
            await page.waitForTimeout(1000);
          },
          { page }
        );
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

      await selectCanonicalMarkerGeneTable(page);
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

  await tryUntil(
    async () => {
      await tooltipLocator.waitFor({ timeout: WAIT_FOR_TIMEOUT_MS });
      const tooltipLocatorVisible = await tooltipLocator.isVisible();
      expect(tooltipLocatorVisible).toBe(true);
    },
    { page }
  );

  // check that tooltip contains text
  const tooltipText = await tooltipLocator.textContent();
  expect(tooltipText).toContain(text);
}

async function selectCanonicalMarkerGeneTable(page: Page) {
  const tableSelector = `[data-testid='${CELL_GUIDE_CARD_CANONICAL_MARKER_GENES_TABLE}']`;

  await tryUntil(
    async () => {
      // set canonical marker genes table as active
      await page
        .getByTestId(CELL_GUIDE_CARD_CANONICAL_MARKER_GENES_TABLE_SELECTOR)
        .click();

      await page
        .locator(tableSelector)
        .waitFor({ timeout: WAIT_FOR_TIMEOUT_MS });
    },
    { page }
  );
}

async function selectComputationalMarkerGeneTable(page: Page) {
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
