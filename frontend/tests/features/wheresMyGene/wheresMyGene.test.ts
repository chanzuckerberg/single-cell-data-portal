import { expect, Page, test, Locator } from "@playwright/test";
import { ROUTES } from "src/common/constants/routes";
import type { RawPrimaryFilterDimensionsResponse } from "src/common/queries/wheresMyGene";
import { FMG_GENE_STRENGTH_THRESHOLD } from "src/views/WheresMyGene/common/constants";
import {
  expandTissue,
  goToPage,
  selectNthOption,
  tryUntil,
  waitForLoadingSpinnerToResolve,
} from "tests/utils/helpers";
import { TEST_URL } from "../../common/constants";
import { TISSUE_DENY_LIST } from "../../fixtures/wheresMyGene/tissueRollup";
import fs from "fs";
import { parse } from "csv-parse/sync";
import AdmZip from "adm-zip";
import {
  WMG_WITH_SEEDED_GENES,
  conditionallyRunTests,
  goToWMG,
  goToWMGWithSeededState,
  waitForHeatmapToRender,
} from "tests/utils/wmgUtils";
import { getCurrentDate } from "tests/utils/downloadUtils";
import { searchAndAddGene } from "tests/utils/geneUtils";

const HOMO_SAPIENS_TERM_ID = "NCBITaxon:9606";

const GENE_LABELS_ID = "[data-testid^=gene-label-]";
const CELL_TYPE_LABELS_ID = "cell-type-name";
const TISSUE_LABELS_ID = "tissue-name";
const ADD_GENE_ID = "add-gene-btn";
const SOURCE_DATA_BUTTON_ID = "source-data-button";
const SOURCE_DATA_LIST_SELECTOR = `[data-testid="source-data-list"]`;
const DOWNLOAD_BUTTON_ID = "download-button";

// FMG test IDs
const ADD_TO_DOT_PLOT_BUTTON_TEST_ID = "add-to-dotplot-fmg-button";
const NO_MARKER_GENES_WARNING_TEST_ID = "no-marker-genes-warning";
const MARKER_SCORES_FMG_TEST_ID = "marker-scores-fmg";

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
const CSV_START_FROM_ROW_NUM = 10; // This is the number of metadata rows + 1
const PNG_CHECKBOX_ID = "png-checkbox";
const CSV_CHECKBOX_ID = "csv-checkbox";
const SVG_CHECKBOX_ID = "svg-checkbox";
const DIALOG_DOWNLOAD_BUTTON_ID = "dialog-download-button";

const MUI_CHIP_ROOT = ".MuiChip-root";

const CELL_TYPE_SANITY_CHECK_NUMBER = 100;

const COMPARE_DROPDOWN_ID = "compare-dropdown";

const EXPORT_OUTPUT_DIR = "playwright-report/";

const FILTERS_PANEL = "filters-panel";

// Error messages
const ERROR_NO_TESTID_OR_LOCATOR = "Either testId or locator must be defined";

const { describe } = test;

describe("Where's My Gene", () => {
  conditionallyRunTests({ forceRun: true });

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
      locator: getDatasetSelectorButton(),
    });
    const numberOfDatasetsBefore = await countLocator(page.getByRole("option"));
    await page.keyboard.press("Escape");

    await clickUntilOptionsShowUp({
      page,
      locator: getDiseaseSelectorButton(),
    });
    const datasetOptions = await page
      .getByRole("option")
      .getByText("acute kidney failure")
      .elementHandles();
    await datasetOptions[0].click();
    await page.keyboard.press("Escape");

    // wait for spinner to disappear, coincides with y-axis update
    await waitForLoadingSpinnerToResolve(page);

    const numberOfTissuesAfter = await countLocator(
      page.getByTestId(TISSUE_LABELS_ID)
    );

    await clickUntilOptionsShowUp({
      page,
      locator: getDatasetSelectorButton(),
    });
    const numberOfDatasetsAfter = await countLocator(page.getByRole("option"));
    await page.keyboard.press("Escape");

    expect(numberOfDatasetsBefore).toBeGreaterThan(numberOfDatasetsAfter);
    expect(numberOfTissuesBefore).toBeGreaterThan(numberOfTissuesAfter);

    function getDiseaseSelector() {
      return page.getByTestId(FILTERS_PANEL).getByTestId("disease-filter");
    }

    function getDiseaseSelectorButton() {
      return getDiseaseSelector().getByRole("button");
    }

    function getDatasetSelector() {
      return page.getByTestId(FILTERS_PANEL).getByTestId("dataset-filter");
    }

    function getDatasetSelectorButton() {
      return getDatasetSelector().getByRole("button");
    }
  });

  test("Source Data", async ({ page }) => {
    await goToWMG(page);

    await clickUntilSidebarShowsUp({ page, testId: SOURCE_DATA_BUTTON_ID });
    await expect(
      page.getByText(
        "After filtering cells with low coverage (less than 500 genes expressed)"
      )
    ).toBeTruthy();

    await tryUntil(
      async () => {
        const numSourceDataListItems = await countLocator(
          page.locator(SOURCE_DATA_LIST_SELECTOR).locator(".MuiListItem-root")
        );
        expect(numSourceDataListItems).toBeGreaterThan(0);

        await page.mouse.click(0, 0);

        function getDatasetSelector() {
          return page.getByTestId(FILTERS_PANEL).getByTestId("dataset-filter");
        }

        const numSelectedDatasetsBefore = await countLocator(
          getDatasetSelector().locator(MUI_CHIP_ROOT)
        );
        expect(numSelectedDatasetsBefore).toBe(0);
        await clickUntilOptionsShowUp({ page, locator: getDatasetSelector() });
        await selectFirstOption(page);
        await clickUntilSidebarShowsUp({ page, testId: SOURCE_DATA_BUTTON_ID });

        const numSourceDataListAfterItems = await countLocator(
          page.locator(SOURCE_DATA_LIST_SELECTOR).locator(".MuiListItem-root")
        );
        expect(numSourceDataListAfterItems).toBeGreaterThan(0);
      },
      { page }
    );
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

      // Check all 3 Compare options work
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

      await clickDropdownOptionByName({
        name: "Sex",
        page,
        testId: COMPARE_DROPDOWN_ID,
      });

      await waitForLoadingSpinnerToResolve(page);

      await tryUntil(
        async () => {
          const afterCellTypeNames = await getCellTypeNames(page);

          expect(
            afterCellTypeNames.find((name) => name.includes("    female"))
          ).toBeTruthy();
        },
        { page }
      );

      await clickDropdownOptionByName({
        name: "Ethnicity",
        page,
        testId: COMPARE_DROPDOWN_ID,
      });

      await waitForLoadingSpinnerToResolve(page);

      await tryUntil(
        async () => {
          const afterCellTypeNames = await getCellTypeNames(page);

          expect(
            afterCellTypeNames.find((name) => name.includes("    multiethnic"))
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

    test("Cell types with no marker genes display warning", async ({
      page,
    }) => {
      await goToWMG(page);

      await expandTissue(page, "heart");

      await getCellTypeFmgButtonAndClick(page, "dendritic cell");

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

  describe("Export CSV ", () => {
    test("Download CSV and validate length - no compare", async ({ page }) => {
      await goToWMG(page);
      const SELECT_N_GENES = 3;

      // Select gene
      await clickUntilOptionsShowUp({ page, testId: ADD_GENE_ID });
      await selectFirstNOptions(SELECT_N_GENES, page);

      await clickUntilDownloadModalShowsUp({
        page,
        testId: DOWNLOAD_BUTTON_ID,
      });

      await page.getByTestId(PNG_CHECKBOX_ID).click();
      await page.getByTestId(CSV_CHECKBOX_ID).click();

      // Start waiting for download before clicking. Note no await.
      const downloadPromise = page.waitForEvent("download");

      await page.getByTestId(DIALOG_DOWNLOAD_BUTTON_ID).click();

      // Wait for download
      const download = await downloadPromise;

      const fileName = download.suggestedFilename();

      expect(fileName).toBe(
        `CELLxGENE_gene_expression_${getCurrentDate()}.csv`
      );

      const filePath = EXPORT_OUTPUT_DIR + fileName;

      await download.saveAs(filePath);

      const csvBuffer = fs.readFileSync(filePath);

      // Parsing will validate rows have consistent column counts per row
      const data = parse(csvBuffer, {
        columns: true,
        from_line: CSV_START_FROM_ROW_NUM, // Start on row with header names, skipping metadata
        skip_empty_lines: true,
      });

      // Validate number of rows
      const cellTypes = await getCellTypeNames(page);
      const tissues = await getTissueNames(page);
      expect(data.length).toBe(
        (cellTypes.length + tissues.length) * SELECT_N_GENES
      );
    });

    test("Download CSV and validate length - compare", async ({ page }) => {
      await goToWMG(page);

      const SELECT_N_GENES = 3;

      // Expand Tissue
      await expandTissue(page, "blood");

      // Select gene
      await clickUntilOptionsShowUp({ page, testId: ADD_GENE_ID });
      await selectFirstNOptions(SELECT_N_GENES, page);

      // Select a compare dimension
      await clickDropdownOptionByName({
        name: "Sex",
        page,
        testId: COMPARE_DROPDOWN_ID,
      });

      await waitForHeatmapToRender(page);

      await clickUntilDownloadModalShowsUp({
        page,
        testId: DOWNLOAD_BUTTON_ID,
      });

      await page.getByTestId(PNG_CHECKBOX_ID).click();
      await page.getByTestId(CSV_CHECKBOX_ID).click();

      // Start waiting for download before clicking. Note no await.
      const downloadPromise = page.waitForEvent("download");

      await page.getByTestId(DIALOG_DOWNLOAD_BUTTON_ID).click();

      // Wait for download
      const download = await downloadPromise;

      const fileName = download.suggestedFilename();

      expect(fileName).toBe(
        `CELLxGENE_gene_expression_${getCurrentDate()}.csv`
      );

      const filePath = EXPORT_OUTPUT_DIR + fileName;

      await download.saveAs(filePath);

      const csvBuffer = fs.readFileSync(filePath);

      // Parsing will validate rows have consistent column counts per row
      const data = parse(csvBuffer, {
        columns: true,
        from_line: CSV_START_FROM_ROW_NUM, // Start on row with header names, skipping metadata
        skip_empty_lines: true,
      });

      // Validate number of rows
      const cellTypes = await getCellTypeNames(page);
      const tissues = await getTissueNames(page);
      expect(data.length).toBe(
        (cellTypes.length + tissues.length) * SELECT_N_GENES
      );
    });
  });

  describe("Multiple file download - one tissue", () => {
    test("Download zip of all outputs (png,svg,csv) for one tissue", async ({
      page,
    }) => {
      await goToWMG(page);

      await expandTissue(page, "blood");

      await clickUntilOptionsShowUp({ page, testId: ADD_GENE_ID });
      await selectFirstOption(page);

      await waitForHeatmapToRender(page);

      await clickUntilDownloadModalShowsUp({
        page,
        testId: DOWNLOAD_BUTTON_ID,
      });

      await page.getByTestId(SVG_CHECKBOX_ID).click();
      await page.getByTestId(CSV_CHECKBOX_ID).click();

      // Start waiting for download before clicking. Note no await.
      const downloadPromise = page.waitForEvent("download");

      await page.getByTestId(DIALOG_DOWNLOAD_BUTTON_ID).click();

      // Wait for download
      const download = await downloadPromise;

      const fileName = download.suggestedFilename();

      expect(fileName).toBe("CELLxGENE_gene_expression.zip");

      const filePath = EXPORT_OUTPUT_DIR + fileName;

      await download.saveAs(filePath);

      const zip = new AdmZip(filePath);

      const zipEntries = zip.getEntries();

      const files = [
        `CELLxGENE_gene_expression_${getCurrentDate()}.csv`,
        `CELLxGENE_gene_expression_${getCurrentDate()}.png`,
        `CELLxGENE_gene_expression_${getCurrentDate()}.svg`,
      ];

      expect(zipEntries.length).toBe(files.length);

      for (const entry of zipEntries) {
        expect(files.includes(entry.name)).toBe(true);
      }
    });
  });

  describe("Multiple file download - two tissues", () => {
    test("Download zip of all outputs (png,svg,csv) for two tissues", async ({
      page,
    }) => {
      await goToWMGWithSeededState(page);

      await waitForHeatmapToRender(page);

      await clickUntilDownloadModalShowsUp({
        page,
        testId: DOWNLOAD_BUTTON_ID,
      });

      await page.getByTestId(SVG_CHECKBOX_ID).click();
      await page.getByTestId(CSV_CHECKBOX_ID).click();

      // Start waiting for download before clicking. Note no await.
      const downloadPromise = page.waitForEvent("download");

      await page.getByTestId(DIALOG_DOWNLOAD_BUTTON_ID).click();

      // Wait for download
      const download = await downloadPromise;

      const fileName = download.suggestedFilename();

      expect(fileName).toBe("CELLxGENE_gene_expression.zip");

      const filePath = "playwright-report/" + fileName;

      await download.saveAs(filePath);

      const zip = new AdmZip(filePath);

      const zipEntries = zip.getEntries();

      const files = [
        `CELLxGENE_gene_expression_${getCurrentDate()}.csv`,
        `CELLxGENE_gene_expression_${getCurrentDate()}.png`,
        `CELLxGENE_gene_expression_${getCurrentDate()}.svg`,
      ];

      expect(zipEntries.length).toBe(files.length);

      for (const entry of zipEntries) {
        expect(files.includes(entry.name)).toBe(true);
      }
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

async function clickUntilDownloadModalShowsUp({
  page,
  testId,
  locator,
}: {
  page: Page;
  testId?: string;
  locator?: Locator;
}) {
  await tryUntil(
    async () => {
      if (testId) {
        await page.getByTestId(testId).click();
      } else if (locator) {
        await locator.click();
      } else {
        throw Error(ERROR_NO_TESTID_OR_LOCATOR);
      }
      await page.locator(".bp4-dialog").elementHandle();
    },
    { page }
  );
}

async function clickUntilSidebarShowsUp({
  page,
  testId,
  locator,
}: {
  page: Page;
  testId?: string;
  locator?: Locator;
}) {
  await tryUntil(
    async () => {
      if (testId) {
        await page.getByTestId(testId).click();
      } else if (locator) {
        await locator.click();
      } else {
        throw Error(ERROR_NO_TESTID_OR_LOCATOR);
      }
      await page.locator(".bp4-drawer-header").elementHandle();
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

async function getTissueNames(page: Page) {
  return getNames({ page, testId: TISSUE_LABELS_ID });
}

// (alec) use this instead of locator.count() to make sure that the element is actually present
// when counting
async function countLocator(locator: Locator) {
  return (await locator.elementHandles()).length;
}
