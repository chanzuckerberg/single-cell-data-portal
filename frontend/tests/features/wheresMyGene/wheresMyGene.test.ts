import { expect, Page, Locator } from "@playwright/test";
import { ROUTES } from "src/common/constants/routes";
import type { RawPrimaryFilterDimensionsResponse } from "src/common/queries/wheresMyGene";
import {
  FMG_GENE_STRENGTH_THRESHOLD,
  FMG_SPECIFICITY_THRESHOLD,
} from "src/views/WheresMyGeneV2/common/constants";
import {
  checkTooltipContent,
  expandTissue,
  goToPage,
  isElementVisible,
  selectNthOption,
  tryUntil,
  waitForLoadingSpinnerToResolve,
} from "tests/utils/helpers";
import { COMPARE_DROPDOWN_ID, TEST_URL } from "../../common/constants";
import { TISSUE_DENY_LIST } from "../../fixtures/wheresMyGene/tissueRollup";

import {
  WMG_WITH_SEEDED_GENES,
  goToWMG,
  goToWMGWithSeededState,
  removeFilteredCellType,
  searchAndAddFilterCellType,
  waitForHeatmapToRender,
} from "tests/utils/wmgUtils";

import {
  CELL_TYPE_SEARCH_PLACEHOLDER_TEXT,
  addGene,
  searchAndAddGene,
} from "tests/utils/geneUtils";
import {
  CELL_TYPE_NAME_LABEL_CLASS_NAME,
  TISSUE_NAME_LABEL_CLASS_NAME,
} from "src/views/WheresMyGeneV2/components/HeatMap/components/YAxisChart/constants";
import { test } from "tests/common/test";
import { MAX_EXPRESSION_LABEL_TEST_ID } from "src/views/WheresMyGeneV2/components/InfoPanel/components/RelativeGeneExpression/constants";
import assert from "assert";
import {
  NO_MARKER_GENES_DESCRIPTION,
  NO_MARKER_GENES_FOR_BLOOD_DESCRIPTION,
  TABLE_HEADER_SPECIFICITY,
  TOO_FEW_CELLS_NO_MARKER_GENES_DESCRIPTION,
} from "src/views/WheresMyGeneV2/components/CellInfoSideBar/constants";
import {
  EFFECT_SIZE,
  MARKER_SCORE_TOOLTIP_CONTENT,
  MARKER_SCORE_TOOLTIP_LINK_TEXT,
  MARKER_SCORE_TOOLTIP_TEST_ID,
  SPECIFICITY_TOOLTIP_CONTENT_FIRST_HALF,
  SPECIFICITY_TOOLTIP_CONTENT_SECOND_HALF,
  SPECIFICITY_TOOLTIP_TEST_ID,
} from "src/common/constants/markerGenes";
import { CELL_GUIDE_CARD_GPT_DESCRIPTION } from "src/views/CellGuide/components/CellGuideCard/components/Description/constants";

const HOMO_SAPIENS_TERM_ID = "NCBITaxon:9606";

const GENE_LABELS_ID = "[data-testid^=gene-label-]";
const CELL_TYPE_LABELS_ID = CELL_TYPE_NAME_LABEL_CLASS_NAME;
const TISSUE_LABELS_ID = TISSUE_NAME_LABEL_CLASS_NAME;
const ADD_GENE_ID = "add-gene-btn";

// FMG test IDs
const ADD_TO_DOT_PLOT_BUTTON_TEST_ID = "add-to-dotplot-fmg-button";
const NO_MARKER_GENES_WARNING_TEST_ID = "no-marker-genes-warning";
const MARKER_SCORES_FMG_TEST_ID = "marker-scores-fmg";
const NO_MARKER_GENES_DESCRIPTION_ID = "no-marker-genes-description";
const EFFECT_SIZE_HEADER_ID = "marker-genes-table-header-score";
const SPECIFICITY_HEADER_ID = "marker-genes-table-specificity";
const SPECIFICITY_TEST_ID = "specificity-fmg";

// gene info test IDs
const GENE_INFO_BUTTON_X_AXIS_TEST_ID = "gene-info-button-x-axis";
const GENE_INFO_HEADER_TEST_ID = "gene-info-header";
const GENE_INFO_GENE_SYNONYMS_TEST_ID = "gene-info-gene-synonyms";
const GENE_INFO_NCBI_LINK_TEST_ID = "gene-info-ncbi-link";
const RIGHT_SIDEBAR_TITLE_TEST_ID = "right-sidebar-title";
const GENE_INFO_TITLE_SPLIT_TEST_ID = "gene-info-title-split";
const GENE_INFO_CLOSE_BUTTON_SPLIT_TEST_ID = "gene-info-close-button-split";
const RIGHT_SIDEBAR_CLOSE_BUTTON_TEST_ID = "right-sidebar-close-button";
const GENE_INFO_BUTTON_CELL_INFO_TEST_ID = "gene-info-button-cell-info";

// Export constants

const MUI_CHIP_ROOT = ".MuiChip-root";

const CELL_TYPE_SANITY_CHECK_NUMBER = 100;

const FILTERS_PANEL = "filters-panel";

// Error messages
const ERROR_NO_TESTID_OR_LOCATOR = "Either testId or locator must be defined";

const { describe } = test;

describe("Where's My Gene", () => {
  test("renders the getting started UI", async ({ page }) => {
    await goToWMG(page);

    await clickUntilOptionsShowUp({ page, testId: ADD_GENE_ID });
    await selectFirstOption(page);

    const filtersPanel = page.getByTestId(FILTERS_PANEL);

    await expect(filtersPanel.getByText("Dataset")).toBeTruthy();
    await expect(filtersPanel.getByText("Disease")).toBeTruthy();
    await expect(
      filtersPanel.getByText("Self-Reported Ethnicity")
    ).toBeTruthy();
    await expect(filtersPanel.getByText("Sex")).toBeTruthy();
    await expect(filtersPanel.getByText("Publication")).toBeTruthy();
    await expect(filtersPanel.getByText("Tissue")).toBeTruthy();

    // Legend
    const Legend = page.getByTestId("legend-wrapper");
    await expect(Legend.getByText("Gene Expression")).toBeTruthy();
    await expect(Legend.getByText("Expressed in Cells (%)")).toBeTruthy();

    // Info Panel
    await expect(filtersPanel.getByText("Methodology")).toBeTruthy();
    await expect(
      filtersPanel.getByText("After filtering cells with low coverage ")
    ).toBeTruthy();
    await expect(filtersPanel.getByText("Source Data")).toBeTruthy();
  });

  test("Filters and Heatmap", async ({ page }) => {
    await goToWMGWithSeededState(page);

    const sexSelector = getSexSelector();
    const selectedSexesBefore = await sexSelector
      .locator(MUI_CHIP_ROOT)
      .elementHandles();

    expect(selectedSexesBefore.length).toBe(0);

    await clickUntilOptionsShowUp({ page, locator: getSexSelectorButton() });

    await selectFirstOption(page);

    const selectedSexesAfter = await sexSelector
      .locator(MUI_CHIP_ROOT)
      .elementHandles();

    expect(selectedSexesAfter.length).toBe(1);

    function getSexSelector() {
      return page.getByTestId(FILTERS_PANEL).getByTestId("sex-filter");
    }

    function getSexSelectorButton() {
      return getSexSelector().getByRole("button");
    }
  });

  test("Primary and secondary filter crossfiltering", async ({ page }) => {
    await goToWMG(page);
    await waitForLoadingSpinnerToResolve(page);

    const numberOfTissuesBefore = await countLocator(
      page.getByTestId(TISSUE_LABELS_ID)
    );

    await clickUntilOptionsShowUp({
      page,
      locator: getDiseaseSelectorButton(page),
    });

    const diseaseOption = await page
      .getByRole("option")
      .getByText("acute kidney failure");

    await diseaseOption.click();
    await page.keyboard.press("Escape");

    // wait for spinner to disappear, coincides with y-axis update
    await waitForLoadingSpinnerToResolve(page);

    const numberOfTissuesAfter = await countLocator(
      page.getByTestId(TISSUE_LABELS_ID)
    );

    expect(numberOfTissuesBefore).toBeGreaterThan(numberOfTissuesAfter);
  });

  test("Hierarchical Clustering", async ({ page }) => {
    await goToWMGWithSeededState(page);

    const beforeGeneNames = await getGeneNames(page);
    const beforeCellTypeNames = await getCellTypeNames(page);

    expect(beforeGeneNames.length).toBe(WMG_WITH_SEEDED_GENES.genes.length);

    // (thuang): Sometimes when API response is slow, we'll not capture all the
    // cell type names, so a sanity check that we expect at least 100 names
    expect(beforeCellTypeNames.length).toBeGreaterThan(
      CELL_TYPE_SANITY_CHECK_NUMBER
    );

    await page.getByTestId("cell-type-sort-dropdown").click();
    await selectNthOption(page, 2);

    await page.getByTestId("gene-sort-dropdown").click();
    await selectNthOption(page, 2);

    const afterGeneNames = await getGeneNames(page);

    const afterCellTypeNames = await getCellTypeNames(page);

    await tryUntil(
      async () => {
        expect(afterGeneNames.length).toBe(beforeGeneNames.length);
        expect(afterCellTypeNames.length).toBe(beforeCellTypeNames.length);

        expect(afterGeneNames).not.toEqual(beforeGeneNames);
        expect(afterCellTypeNames).not.toEqual(beforeCellTypeNames);
      },
      { page }
    );
  });

  describe("tissue rollup", () => {
    test("does NOT show tissues on the deny list", async ({ page }) => {
      const [response] = await Promise.all([
        page.waitForResponse("**/primary_filter_dimensions"),
        goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page),
      ]);

      const primaryFilterDimensions =
        (await response.json()) as RawPrimaryFilterDimensionsResponse;

      const humanTissues =
        primaryFilterDimensions.tissue_terms[HOMO_SAPIENS_TERM_ID];

      const tissueIds = new Set(
        humanTissues.map((tissue) => Object.keys(tissue)[0])
      );

      const hasDeniedTissue = TISSUE_DENY_LIST.some((deniedTissueId) =>
        tissueIds.has(deniedTissueId)
      );

      expect(hasDeniedTissue).toBe(false);
    });
  });

  describe("Compare", () => {
    test("Display stratified labels in y-axis as expected", async ({
      page,
    }) => {
      await goToWMGWithSeededState(page);
      await waitForHeatmapToRender(page);

      await tryUntil(
        async () => {
          const beforeCellTypeNames = await getCellTypeNames(page);

          // (thuang): Sometimes when API response is slow, we'll not capture all the
          // cell type names, so a sanity check that we expect at least 100 names
          expect(beforeCellTypeNames.length).toBeGreaterThan(
            CELL_TYPE_SANITY_CHECK_NUMBER
          );

          // beforeCellTypeNames array does not contain "normal"

          /**
           * (thuang): Make sure the default y axis is not stratified by checking
           * that there are no 4 spaces in the cell type names for indentation
           */
          expect(
            beforeCellTypeNames.find((name) => name.includes("    "))
          ).toBeFalsy();
        },
        { page }
      );

      // Check disease compare option works
      // (alec) Checking all compare options is too slow, so just check one
      // the risk of the other compare options deviating from the behavior of the
      // disease option is low
      // Wait for loading spinner to disappear before checking (compare data is slow to load)
      await clickDropdownOptionByName({
        name: "Disease",
        page,
        testId: COMPARE_DROPDOWN_ID,
      });

      await waitForLoadingSpinnerToResolve(page);

      await tryUntil(
        async () => {
          const afterCellTypeNames = await getCellTypeNames(page);
          expect(
            afterCellTypeNames.find((name) => name.includes("    normal"))
          ).toBeTruthy();
        },
        { page }
      );

      // Check selecting None removes stratification
      await clickDropdownOptionByName({
        name: "None",
        page,
        testId: COMPARE_DROPDOWN_ID,
      });

      await tryUntil(
        async () => {
          expect(
            (await getCellTypeNames(page)).find((name) => name.includes("    "))
          ).toBeFalsy();
        },
        { page }
      );
    });
  });

  describe("Find Marker Genes", () => {
    test("Marker gene panel opens and can add genes to dotplot", async ({
      page,
    }) => {
      await goToWMG(page);

      await expandTissue(page, "lung");

      await getCellTypeFmgButtonAndClick(page, "memory B cell");

      await getButtonAndClick(page, ADD_TO_DOT_PLOT_BUTTON_TEST_ID);

      await waitForHeatmapToRender(page);

      const geneLabels = await getGeneNames(page);
      const numGenes = geneLabels.length;
      expect(numGenes).toBeGreaterThan(0);
    });

    test("Marker gene panel displays cell type description", async ({
      page,
    }) => {
      await goToWMG(page);

      await expandTissue(page, "lung");

      await getCellTypeFmgButtonAndClick(page, "memory B cell");

      await waitForElement(page, CELL_GUIDE_CARD_GPT_DESCRIPTION);
    });

    // need to find a tissue, cell type with no marker genes
    test.skip("Cell types with no marker genes display warning", async ({
      page,
    }) => {
      await goToWMG(page);

      await expandTissue(page, "brain");

      await getCellTypeFmgButtonAndClick(page, "smooth muscle cell");

      await waitForElement(page, NO_MARKER_GENES_WARNING_TEST_ID);
    });

    test(`Marker scores are always greater than or equal to ${FMG_GENE_STRENGTH_THRESHOLD}`, async ({
      page,
    }) => {
      await goToWMG(page);

      await expandTissue(page, "lung");

      await getCellTypeFmgButtonAndClick(page, "club cell");

      const markerScores = await page
        .getByTestId(MARKER_SCORES_FMG_TEST_ID)
        .allTextContents();

      for (const markerScore of markerScores) {
        expect(parseFloat(markerScore)).toBeGreaterThanOrEqual(
          FMG_GENE_STRENGTH_THRESHOLD
        );
      }
    });

    // Note: This test could fail if we add more adipose tissue naive B cells to our corpus.
    // If this happens, just mark the test as skipped and ping @joyceyan
    test(`Should verify cell types with < 25 cells have no marker genes`, async ({
      page,
    }) => {
      await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

      // Expand blood tissue
      await expandTissue(page, "adipose-tissue");

      // Click naive B cell info icon
      const naiveBCell = page.getByTestId(
        "cell-type-info-button-adipose tissue-naive B cell"
      );
      await naiveBCell.scrollIntoViewIfNeeded();
      await naiveBCell.click();

      // Verify copy is what we expect
      const noMarkerGenesDescription = (await page
        .getByTestId(NO_MARKER_GENES_DESCRIPTION_ID)
        .textContent()) as string;
      assert.strictEqual(
        noMarkerGenesDescription.trim(),
        TOO_FEW_CELLS_NO_MARKER_GENES_DESCRIPTION
      );
    });

    // Note: This test could fail if we find more marker genes for embryo hematopoietic cells.
    // If this happens, just mark the test as skipped and ping @joyceyan
    test(`Should verify copy for cell types with no marker genes`, async ({
      page,
    }) => {
      await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

      // Expand embryo tissue
      await expandTissue(page, "embryo");

      // Click yolk sac somatic cell info icon
      const hematopoieticCell = page.getByTestId(
        "cell-type-info-button-embryo-hematopoietic cell"
      );
      await hematopoieticCell.scrollIntoViewIfNeeded();
      await hematopoieticCell.click();

      // Verify copy is what we expect
      const noMarkerGenesDescription = (await page
        .getByTestId(NO_MARKER_GENES_DESCRIPTION_ID)
        .textContent()) as string;
      assert.strictEqual(
        noMarkerGenesDescription.trim(),
        NO_MARKER_GENES_DESCRIPTION
      );
    });

    test(`Should verify blood cells have no marker genes`, async ({ page }) => {
      await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

      // Expand blood tissue
      await expandTissue(page, "blood");

      // Click stem cell info icon
      await page.getByTestId("cell-type-info-button-blood-stem cell").click();

      // Verify copy is what we expect
      const noMarkerGenesDescription = (await page
        .getByTestId(NO_MARKER_GENES_DESCRIPTION_ID)
        .textContent()) as string;
      assert.strictEqual(
        noMarkerGenesDescription.trim(),
        NO_MARKER_GENES_FOR_BLOOD_DESCRIPTION
      );
    });

    test(`Should verify effect size and specificity column`, async ({
      page,
    }) => {
      await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

      // Expand adipose tissue
      await expandTissue(page, "adipose-tissue");

      // Click into a cell type that has marker genes
      await page
        .getByTestId("cell-type-info-button-adipose tissue-contractile cell")
        .click();

      // Verify effect size header and tooltip
      const effectSizeHeader = (await page
        .getByTestId(EFFECT_SIZE_HEADER_ID)
        .textContent()) as string;
      assert.strictEqual(effectSizeHeader.trim(), EFFECT_SIZE);

      await isElementVisible(page, MARKER_SCORE_TOOLTIP_TEST_ID);
      await page.getByTestId(MARKER_SCORE_TOOLTIP_TEST_ID).hover();
      await checkTooltipContent(page, MARKER_SCORE_TOOLTIP_CONTENT);
      await checkTooltipContent(page, MARKER_SCORE_TOOLTIP_LINK_TEXT);

      // Verify specificity header and tooltip
      const specificityHeader = (await page
        .getByTestId(SPECIFICITY_HEADER_ID)
        .textContent()) as string;
      assert.strictEqual(specificityHeader.trim(), TABLE_HEADER_SPECIFICITY);

      await isElementVisible(page, SPECIFICITY_TOOLTIP_TEST_ID);
      await page.getByTestId(SPECIFICITY_TOOLTIP_TEST_ID).hover();
      await checkTooltipContent(
        page,
        SPECIFICITY_TOOLTIP_CONTENT_FIRST_HALF +
          " adipose tissue " +
          SPECIFICITY_TOOLTIP_CONTENT_SECOND_HALF
      );

      // Verify specificity values are valid
      const specificityScores = await page
        .getByTestId(SPECIFICITY_TEST_ID)
        .allTextContents();

      for (const specificityScore of specificityScores) {
        expect(parseFloat(specificityScore)).toBeGreaterThanOrEqual(
          FMG_SPECIFICITY_THRESHOLD
        );
      }
    });
  });

  describe("Gene info", () => {
    const TEST_GENE = "DMP1";

    test("Display gene info panel in sidebar", async ({ page }) => {
      await searchAndAddGene(page, TEST_GENE);

      await waitForHeatmapToRender(page);

      // hover over gene label
      await page.getByTestId(`gene-name-${TEST_GENE}`).hover();

      await getFirstButtonAndClick(page, GENE_INFO_BUTTON_X_AXIS_TEST_ID);

      await waitForElement(page, GENE_INFO_HEADER_TEST_ID);
      await waitForElement(page, GENE_INFO_GENE_SYNONYMS_TEST_ID);
      await waitForElement(page, GENE_INFO_NCBI_LINK_TEST_ID);
      await waitForElement(page, RIGHT_SIDEBAR_TITLE_TEST_ID);

      await getButtonAndClick(page, RIGHT_SIDEBAR_CLOSE_BUTTON_TEST_ID);

      await waitForElementToBeRemoved(page, RIGHT_SIDEBAR_TITLE_TEST_ID);
    });

    test("Display gene info bottom drawer in cell info sidebar", async ({
      page,
    }) => {
      await searchAndAddGene(page, TEST_GENE);

      await expandTissue(page, "lung");

      await waitForHeatmapToRender(page);

      // hover over gene label
      await page.getByTestId(`gene-name-${TEST_GENE}`).hover();

      await waitForHeatmapToRender(page);

      await getCellTypeFmgButtonAndClick(page, "muscle cell");

      await getFirstButtonAndClick(page, GENE_INFO_BUTTON_CELL_INFO_TEST_ID);

      await waitForElement(page, GENE_INFO_HEADER_TEST_ID);
      await waitForElement(page, GENE_INFO_GENE_SYNONYMS_TEST_ID);
      await waitForElement(page, GENE_INFO_NCBI_LINK_TEST_ID);
      await waitForElement(page, GENE_INFO_TITLE_SPLIT_TEST_ID);

      await getButtonAndClick(page, GENE_INFO_CLOSE_BUTTON_SPLIT_TEST_ID);

      await waitForElementToBeRemoved(page, GENE_INFO_TITLE_SPLIT_TEST_ID);

      await getButtonAndClick(page, RIGHT_SIDEBAR_CLOSE_BUTTON_TEST_ID);

      await waitForElementToBeRemoved(page, RIGHT_SIDEBAR_TITLE_TEST_ID);

      // hover over gene label
      await page.getByTestId(`gene-name-${TEST_GENE}`).hover();

      await getFirstButtonAndClick(page, GENE_INFO_BUTTON_X_AXIS_TEST_ID);

      await waitForElement(page, GENE_INFO_HEADER_TEST_ID);
      await waitForElement(page, RIGHT_SIDEBAR_TITLE_TEST_ID);
    });
  });

  describe("Newsletter", () => {
    const NEWSLETTER_MODAL_CONTENT = "newsletter-modal-content";
    const NEWSLETTER_MODAL_OPEN_BUTTON = "newsletter-modal-open-button";
    const NEWSLETTER_MODAL_CLOSE_BUTTON = "newsletter-modal-close-button";
    const NEWSLETTER_SUBSCRIBE_BUTTON = "newsletter-subscribe-button";
    const NEWSLETTER_EMAIL_INPUT = "newsletter-email-input";
    const NEWSLETTER_VALIDATION_ERROR_MESSAGE =
      "newsletter-validation-error-message";
    const FAILED_EMAIL_VALIDATION_STRING =
      "Please provide a valid email address.";

    test("Newsletter Modal - Open/Close", async ({ page }) => {
      await goToWMG(page);

      // Open modal
      await getButtonAndClick(page, NEWSLETTER_MODAL_OPEN_BUTTON);

      await waitForElement(page, NEWSLETTER_MODAL_CONTENT);

      const modalContent = page.getByTestId(NEWSLETTER_MODAL_CONTENT);

      // modal content
      expect(modalContent.getByText("Join Our Newsletter")).toBeTruthy();
      expect(
        modalContent.getByText(
          "Get a quarterly email with the latest CELLxGENE features and data."
        )
      ).toBeTruthy();
      expect(modalContent.getByText("Enter email address")).toBeTruthy();
      expect(modalContent.getByText("Subscribe")).toBeTruthy();
      expect(modalContent.getByText("Unsubscribe at any time.")).toBeTruthy();

      // Close modal
      await getButtonAndClick(page, NEWSLETTER_MODAL_CLOSE_BUTTON);

      await waitForElementToBeRemoved(page, NEWSLETTER_MODAL_CONTENT);
    });

    test("Newsletter Modal - Validate Email", async ({ page }) => {
      await goToWMG(page);

      // Open modal
      await getButtonAndClick(page, NEWSLETTER_MODAL_OPEN_BUTTON);

      await waitForElement(page, NEWSLETTER_MODAL_CONTENT);

      const emailInput = page.getByTestId(NEWSLETTER_EMAIL_INPUT);
      const subscribeButton = page.getByTestId(NEWSLETTER_SUBSCRIBE_BUTTON);
      const validationMessage = page.getByTestId(
        NEWSLETTER_VALIDATION_ERROR_MESSAGE
      );

      // No input
      expect(subscribeButton.isDisabled());

      // Bad email 1
      emailInput.fill("test");
      expect(subscribeButton.isEnabled());
      await subscribeButton.click();
      expect(
        validationMessage.getByText(FAILED_EMAIL_VALIDATION_STRING)
      ).toBeTruthy();

      emailInput.fill("");
      expect(subscribeButton.isDisabled());

      // Bad email 2
      emailInput.fill("test@test");
      expect(subscribeButton.isEnabled());
      await subscribeButton.click();
      expect(
        validationMessage.getByText(FAILED_EMAIL_VALIDATION_STRING)
      ).toBeTruthy();

      emailInput.fill("");
      expect(subscribeButton.isDisabled());

      // Bad email 3
      emailInput.fill("test@test.chanzuckerberg");
      expect(subscribeButton.isEnabled());
      await subscribeButton.click();
      expect(
        validationMessage.getByText(FAILED_EMAIL_VALIDATION_STRING)
      ).toBeTruthy();
    });
  });

  describe("Clear All Genes Button", () => {
    const CLEAR_GENES_BUTTON_ID = "clear-genes-button";

    test("Clear three genes", async ({ page }) => {
      await goToWMGWithSeededState(page);

      // Genes before clear
      const beforeGeneNames = await getGeneNames(page);
      expect(beforeGeneNames.length).toBe(WMG_WITH_SEEDED_GENES.genes.length);

      // Click clear all button
      await page.getByTestId(CLEAR_GENES_BUTTON_ID).click();

      // Count genes after clear
      const afterGeneNames = await getGeneNames(page);

      await tryUntil(
        async () => {
          expect(afterGeneNames.length).toBe(0);
          expect(afterGeneNames).not.toEqual(beforeGeneNames);
        },
        { page }
      );
    });
  });

  /**
   * https://github.com/chanzuckerberg/single-cell-data-portal/issues/5715
   */
  test("Bug fix #5715: Adding a gene collapses open tissues", async ({
    page,
  }) => {
    await goToWMG(page);
    await expandTissue(page, "lung");

    await waitForLoadingSpinnerToResolve(page);

    await tryUntil(
      async () => {
        expect((await getCellTypeNames(page)).length).toBeGreaterThan(0);
      },
      { page }
    );

    const beforeAddGeneCellTypeCount = (await getCellTypeNames(page)).length;

    await Promise.all([
      waitForLoadingSpinnerToResolve(page),
      addGene(page, "MALAT1"),
    ]);

    await tryUntil(
      async () => {
        const afterAddGeneCellTypeCount = (await getCellTypeNames(page)).length;

        expect(afterAddGeneCellTypeCount).toBe(beforeAddGeneCellTypeCount);
      },
      { page }
    );
  });
  describe("Cell Type Filtering", () => {
    test("Cell Types don't contain UBERON terms", async ({ page }) => {
      await goToWMG(page);
      await waitForLoadingSpinnerToResolve(page);
      await page
        .getByPlaceholder(CELL_TYPE_SEARCH_PLACEHOLDER_TEXT)
        .fill("UBERON");
      await page.keyboard.type("UBERON");
    });
    test("Filter to multiple cell types and then clear", async ({ page }) => {
      const CELL_TYPE_NAMES = ["B cell", "T cell", "PP cell"];

      await goToWMG(page);
      await waitForLoadingSpinnerToResolve(page);

      for (const cellTypeName of CELL_TYPE_NAMES) {
        await searchAndAddFilterCellType(page, cellTypeName);
      }

      for (const cellTypeName of CELL_TYPE_NAMES) {
        await removeFilteredCellType(page, cellTypeName);
      }
    });
  });
  describe("Legend dynamic scaling", () => {
    test("Filter to multiple cell types and then clear", async ({ page }) => {
      const CELL_TYPE_NAMES = ["plasma cell"];

      await goToWMG(page);
      await waitForLoadingSpinnerToResolve(page);
      await page
        .getByTestId("newsletter-modal-banner-wrapper")
        .getByLabel("Close")
        .click();
      await clickUntilOptionsShowUp({ page, testId: ADD_GENE_ID });

      await page.keyboard.type("JCHAIN");
      await page.keyboard.press("ArrowDown");
      await page.keyboard.press("Enter");

      await clickDropdownOptionByName({
        page,
        testId: "color-scale-dropdown",
        name: "Unscaled",
      });

      const textContentBefore = await page
        .getByTestId(MAX_EXPRESSION_LABEL_TEST_ID)
        .textContent();
      for (const cellTypeName of CELL_TYPE_NAMES) {
        await searchAndAddFilterCellType(page, cellTypeName);
      }
      const textContentAfter = await page
        .getByTestId(MAX_EXPRESSION_LABEL_TEST_ID)
        .textContent();
      for (const cellTypeName of CELL_TYPE_NAMES) {
        await removeFilteredCellType(page, cellTypeName);
      }
      const textContentBefore2 = await page
        .getByTestId(MAX_EXPRESSION_LABEL_TEST_ID)
        .textContent();

      expect(textContentBefore).not.toEqual(textContentAfter);
      expect(textContentBefore).toEqual(textContentBefore2);
    });
  });
});

async function getNames({
  page,
  testId,
  locator,
}: {
  page: Page;
  testId?: string;
  locator?: Locator;
}): Promise<string[]> {
  let labelsLocator: Locator;
  if (testId) {
    labelsLocator = page.getByTestId(testId);
  } else if (locator) {
    labelsLocator = locator;
  } else {
    throw Error(ERROR_NO_TESTID_OR_LOCATOR);
  }
  await tryUntil(
    async () => {
      const names = await labelsLocator.allTextContents();
      if (names.length) {
        expect(typeof names[0]).toBe("string");
      }
    },
    { page }
  );

  return labelsLocator.allTextContents();
}
function getDiseaseSelector(page: Page) {
  return page.getByTestId(FILTERS_PANEL).getByTestId("disease-filter");
}

function getDiseaseSelectorButton(page: Page) {
  return getDiseaseSelector(page).getByRole("button");
}

async function clickUntilOptionsShowUp({
  page,
  testId,
  locator,
}: {
  page: Page;
  testId?: string;
  locator?: Locator;
}) {
  // either testId or locator must be defined, not both
  // locator is used when the element cannot be found using just the test Id from the page
  await tryUntil(
    async () => {
      if (testId) {
        await page.getByTestId(testId).click();
      } else if (locator) {
        await locator.click();
      } else {
        throw Error(ERROR_NO_TESTID_OR_LOCATOR);
      }
      await page.getByRole("tooltip").getByRole("option").elementHandles();
    },
    { page }
  );
}

// (thuang): This only works when a dropdown is open
export async function selectFirstOption(page: Page) {
  await selectFirstNOptions(1, page);
}

async function selectFirstNOptions(count: number, page: Page) {
  for (let i = 0; i < count; i++) {
    await page.keyboard.press("ArrowDown");
    await page.keyboard.press("Enter");
  }

  await page.keyboard.press("Escape");
}

async function waitForElement(page: Page, testId: string) {
  await tryUntil(
    async () => {
      await expect(page.getByTestId(testId)).not.toHaveCount(0);
    },
    { page }
  );
}

async function waitForElementToBeRemoved(page: Page, testId: string) {
  await tryUntil(
    async () => {
      await expect(page.getByTestId(testId)).toHaveCount(0);
    },
    { page }
  );
}

async function getButtonAndClick(page: Page, testID: string) {
  await tryUntil(
    async () => {
      await page.getByTestId(testID).click();
    },
    { page }
  );
}

// for when there are multiple buttons with the same testID
async function getFirstButtonAndClick(page: Page, testID: string) {
  await tryUntil(
    async () => {
      const buttons = await page.getByTestId(testID).elementHandles();
      await buttons[0].click();
    },
    { page }
  );
}

async function getCellTypeFmgButtonAndClick(page: Page, cellType: string) {
  await waitForElement(page, CELL_TYPE_LABELS_ID);

  await tryUntil(
    async () => {
      await page
        .getByRole("img", {
          name: "display marker genes for " + cellType,
        })
        .click();
    },
    { page }
  );
}

async function clickDropdownOptionByName({
  page,
  testId,
  name,
}: {
  page: Page;
  testId: string;
  name: string;
}) {
  await page.getByTestId(testId).click();
  await page.getByRole("option").filter({ hasText: name }).click();
  await page.keyboard.press("Escape");
}

async function getGeneNames(page: Page) {
  return getNames({ page, locator: page.locator(GENE_LABELS_ID) });
}

async function getCellTypeNames(page: Page) {
  return getNames({ page, testId: CELL_TYPE_LABELS_ID });
}

// (alec) use this instead of locator.count() to make sure that the element is actually present
// when counting
async function countLocator(locator: Locator) {
  return (await locator.elementHandles()).length;
}
