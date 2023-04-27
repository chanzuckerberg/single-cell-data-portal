import { expect, test } from "@playwright/test";
import { ROUTES } from "src/common/constants/routes";
import type { RawPrimaryFilterDimensionsResponse } from "src/common/queries/wheresMyGene";
import { FMG_GENE_STRENGTH_THRESHOLD } from "src/views/WheresMyGene/common/constants";
import {
  clickDropdownOptionByName,
  clickUntilDownloadModalShowsUp,
  clickUntilOptionsShowUp,
  clickUntilSidebarShowsUp,
  countLocator,
  getButtonAndClick,
  getCellTypeFmgButtonAndClick,
  getCellTypeNames,
  getFirstButtonAndClick,
  getGeneNames,
  goToPage,
  isDevStagingProd,
  selectFirstNOptions,
  selectFirstOption,
  selectNthOption,
  tryUntil,
  waitForElement,
  waitForElementToBeRemoved,
  waitForHeatmapToRender,
} from "tests/utils/helpers";
import { TEST_URL } from "../../common/constants";
import { TISSUE_DENY_LIST } from "../../fixtures/wheresMyGene/tissueRollup";
import fs from "fs";
import { parse } from "csv-parse/sync";
import AdmZip from "adm-zip";
import {
  ADD_TISSUE_ID,
  ADD_GENE_ID,
  CELL_TYPE_LABELS_ID,
  SOURCE_DATA_BUTTON_ID,
  SOURCE_DATA_LIST_SELECTOR,
  GENE_DELETE_BUTTON,
  HOMO_SAPIENS_TERM_ID,
  DOWNLOAD_BUTTON_ID,
} from "tests/utils/constants";

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
const CSV_START_FROM_ROW_NUM = 9; // This is the number of metadata rows + 1
const PNG_CHECKBOX_ID = "png-checkbox";
const CSV_CHECKBOX_ID = "csv-checkbox";
const SVG_CHECKBOX_ID = "svg-checkbox";
const DIALOG_DOWNLOAD_BUTTON_ID = "dialog-download-button";

const MUI_CHIP_ROOT = ".MuiChip-root";

const CELL_TYPE_SANITY_CHECK_NUMBER = 100;

const COMPARE_DROPDOWN_ID = "compare-dropdown";

const EXPORT_OUTPUT_DIR = "playwright-report/";

const FILTERS_PANEL = "filters-panel";

const { describe, skip } = test;

describe("Where's My Gene", () => {
  skip(!isDevStagingProd, "WMG BE API does not work locally or in rdev");

  test("renders the getting started UI", async ({ page }) => {
    await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

    // Getting Started section
    expect(page.getByText("STEP 1")).toBeTruthy();
    expect(page.getByText("Add Tissues")).toBeTruthy();

    expect(page.getByText("STEP 2")).toBeTruthy();
    expect(page.getByText("Add Genes")).toBeTruthy();

    expect(page.getByText("STEP 3")).toBeTruthy();
    expect(page.getByText("Explore Gene Expression")).toBeTruthy();

    clickUntilOptionsShowUp({ page, testId: ADD_TISSUE_ID });
    selectFirstOption(page);

    clickUntilOptionsShowUp({ page, testId: ADD_GENE_ID });
    selectFirstOption(page);

    const filtersPanel = page.getByTestId(FILTERS_PANEL);

    expect(filtersPanel.getByText("Dataset")).toBeTruthy();
    expect(filtersPanel.getByText("Disease")).toBeTruthy();
    expect(filtersPanel.getByText("Self-Reported Ethnicity")).toBeTruthy();
    expect(filtersPanel.getByText("Sex")).toBeTruthy();

    // Legend
    const Legend = page.getByTestId("legend-wrapper");
    expect(Legend.getByText("Gene Expression")).toBeTruthy();
    expect(Legend.getByText("Expressed in Cells (%)")).toBeTruthy();

    // Info Panel
    expect(filtersPanel.getByText("Methodology")).toBeTruthy();
    expect(
      filtersPanel.getByText("After filtering cells with low coverage ")
    ).toBeTruthy();
    expect(filtersPanel.getByText("Source Data")).toBeTruthy();
  });

  test("Filters and Heatmap", async ({ page }) => {
    await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

    clickUntilOptionsShowUp({ page, testId: ADD_TISSUE_ID });
    selectFirstOption(page);

    clickUntilOptionsShowUp({ page, testId: ADD_GENE_ID });
    selectFirstOption(page);

    waitForHeatmapToRender(page);

    const sexSelector = getSexSelector();
    const selectedSexesBefore = await sexSelector
      .locator(MUI_CHIP_ROOT)
      .elementHandles();

    expect(selectedSexesBefore.length).toBe(0);

    clickUntilOptionsShowUp({ page, locator: getSexSelectorButton() });

    selectFirstOption(page);

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
    await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

    clickUntilOptionsShowUp({ page, testId: ADD_TISSUE_ID });
    const numberOfTissuesBefore = countLocator(page.getByRole("option"));
    await page.keyboard.press("Escape");

    clickUntilOptionsShowUp({
      page,
      locator: getDiseaseSelectorButton(),
    });
    const numberOfDiseasesBefore = countLocator(page.getByRole("option"));
    await page.keyboard.press("Escape");

    clickUntilOptionsShowUp({
      page,
      locator: getDatasetSelectorButton(),
    });
    const datasetOptions = await page.getByRole("option").elementHandles();
    await datasetOptions[0].click();
    await page.keyboard.press("Escape");

    clickUntilOptionsShowUp({ page, testId: ADD_TISSUE_ID });
    const numberOfTissuesAfter = await countLocator(page.getByRole("option"));
    await page.keyboard.press("Escape");

    clickUntilOptionsShowUp({
      page,
      locator: getDiseaseSelectorButton(),
    });
    const numberOfDiseasesAfter = await countLocator(page.getByRole("option"));
    await page.keyboard.press("Escape");

    expect(numberOfDiseasesBefore).toBeGreaterThan(numberOfDiseasesAfter);
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

  test("Selected tissue and no genes displays cell types", async ({ page }) => {
    await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

    clickUntilOptionsShowUp({ page, testId: ADD_TISSUE_ID });
    selectFirstOption(page);

    waitForElement(page, CELL_TYPE_LABELS_ID);
    const numCellTypes = countLocator(page.getByTestId(CELL_TYPE_LABELS_ID));
    expect(numCellTypes).toBeGreaterThan(0);
  });

  test("Source Data", async ({ page }) => {
    await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

    clickUntilOptionsShowUp({ page, testId: ADD_TISSUE_ID });
    selectFirstOption(page);

    clickUntilOptionsShowUp({ page, testId: ADD_GENE_ID });
    selectFirstOption(page);

    waitForHeatmapToRender(page);
    clickUntilSidebarShowsUp({ page, testId: SOURCE_DATA_BUTTON_ID });
    expect(
      page.getByText(
        "After filtering cells with low coverage (less than 500 genes expressed)"
      )
    ).toBeTruthy();

    await tryUntil(
      async () => {
        const numSourceDataListItems = countLocator(
          page.locator(SOURCE_DATA_LIST_SELECTOR).locator(".MuiListItem-root")
        );
        expect(numSourceDataListItems).toBeGreaterThan(0);

        await page.mouse.click(0, 0);

        function getDatasetSelector() {
          return page.getByTestId(FILTERS_PANEL).getByTestId("dataset-filter");
        }

        const numSelectedDatasetsBefore = countLocator(
          getDatasetSelector().locator(MUI_CHIP_ROOT)
        );
        expect(numSelectedDatasetsBefore).toBe(0);
        clickUntilOptionsShowUp({ page, locator: getDatasetSelector() });
        selectFirstOption(page);
        clickUntilSidebarShowsUp({ page, testId: SOURCE_DATA_BUTTON_ID });

        const numSourceDataListAfterItems = countLocator(
          page.locator(SOURCE_DATA_LIST_SELECTOR).locator(".MuiListItem-root")
        );
        expect(numSourceDataListAfterItems).toBeGreaterThan(0);
      },
      { page }
    );
  });

  test("Hierarchical Clustering", async ({ page }) => {
    await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

    const TISSUE_COUNT = 1;
    const GENE_COUNT = 3;

    clickUntilOptionsShowUp({ page, testId: ADD_TISSUE_ID });
    selectFirstNOptions(TISSUE_COUNT, page);

    clickUntilOptionsShowUp({ page, testId: ADD_GENE_ID });
    selectFirstNOptions(GENE_COUNT, page);

    const beforeGeneNames = getGeneNames(page);
    const beforeCellTypeNames = getCellTypeNames(page);

    expect((await beforeGeneNames).length).toBe(GENE_COUNT);

    // (thuang): Sometimes when API response is slow, we'll not capture all the
    // cell type names, so a sanity check that we expect at least 100 names
    expect((await beforeCellTypeNames).length).toBeGreaterThan(
      CELL_TYPE_SANITY_CHECK_NUMBER
    );

    await page.getByTestId("cell-type-sort-dropdown").click();
    await selectNthOption(page, 2);

    await page.getByTestId("gene-sort-dropdown").click();
    await selectNthOption(page, 2);

    const afterGeneNames = getGeneNames(page);

    const afterCellTypeNames = getCellTypeNames(page);

    await tryUntil(
      async () => {
        expect((await afterGeneNames).length).toBe(
          (await beforeGeneNames).length
        );
        expect((await afterCellTypeNames).length).toBe(
          (await beforeCellTypeNames).length
        );

        expect(afterGeneNames).not.toEqual(beforeGeneNames);
        expect(afterCellTypeNames).not.toEqual(beforeCellTypeNames);
      },
      { page }
    );
  });

  test("delete genes", async ({ page }) => {
    await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

    clickUntilOptionsShowUp({ page, testId: ADD_TISSUE_ID });
    selectFirstNOptions(1, page);

    clickUntilOptionsShowUp({ page, testId: ADD_GENE_ID });
    selectFirstNOptions(3, page);

    waitForHeatmapToRender(page);

    const beforeGeneNames = getGeneNames(page);

    await tryUntil(
      async () => {
        await page.keyboard.press("Backspace");

        // Testing single gene delete
        await page.hover(".gene-label-container");
        getFirstButtonAndClick(page, GENE_DELETE_BUTTON);

        const afterGeneNames = getGeneNames(page);

        expect((await afterGeneNames).length).toBe(
          (await beforeGeneNames).length - 1
        );
        expect(afterGeneNames).not.toEqual(beforeGeneNames);
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
      await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

      clickUntilOptionsShowUp({ page, testId: ADD_TISSUE_ID });
      selectFirstNOptions(1, page);

      clickUntilOptionsShowUp({ page, testId: ADD_GENE_ID });
      selectFirstNOptions(3, page);

      waitForHeatmapToRender(page);

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
           * that there are no 2 spaces in the cell type names for indentation
           */
          expect(
            beforeCellTypeNames.find((name) => name.includes("  "))
          ).toBeFalsy();
        },
        { page }
      );

      // Check all 3 Compare options work
      clickDropdownOptionByName({
        name: "Disease",
        page,
        testId: COMPARE_DROPDOWN_ID,
      });

      await tryUntil(
        async () => {
          const afterCellTypeNames = await getCellTypeNames(page);

          expect(
            afterCellTypeNames.find((name) => name.includes("  normal"))
          ).toBeTruthy();
        },
        { page }
      );

      clickDropdownOptionByName({
        name: "Sex",
        page,
        testId: COMPARE_DROPDOWN_ID,
      });

      await tryUntil(
        async () => {
          const afterCellTypeNames = await getCellTypeNames(page);

          expect(
            afterCellTypeNames.find((name) => name.includes("  female"))
          ).toBeTruthy();
        },
        { page }
      );

      clickDropdownOptionByName({
        name: "Ethnicity",
        page,
        testId: COMPARE_DROPDOWN_ID,
      });

      await tryUntil(
        async () => {
          const afterCellTypeNames = await getCellTypeNames(page);

          expect(
            afterCellTypeNames.find((name) => name.includes("  multiethnic"))
          ).toBeTruthy();
        },
        { page }
      );

      // Check selecting None removes stratification
      clickDropdownOptionByName({
        name: "None",
        page,
        testId: COMPARE_DROPDOWN_ID,
      });

      await tryUntil(
        async () => {
          expect(
            (await getCellTypeNames(page)).find((name) => name.includes("  "))
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
      await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

      clickDropdownOptionByName({
        page,
        testId: ADD_TISSUE_ID,
        name: "lung",
      });

      getCellTypeFmgButtonAndClick(page, "memory B cell");

      getButtonAndClick(page, ADD_TO_DOT_PLOT_BUTTON_TEST_ID);

      waitForHeatmapToRender(page);

      const geneLabels = await getGeneNames(page);
      const numGenes = geneLabels.length;
      expect(numGenes).toBeGreaterThan(0);
    });

    test("Cell types with no marker genes display warning", async ({
      page,
    }) => {
      await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

      clickDropdownOptionByName({
        page,
        testId: ADD_TISSUE_ID,
        name: "heart",
      });

      await getCellTypeFmgButtonAndClick(page, "dendritic cell");

      waitForElement(page, NO_MARKER_GENES_WARNING_TEST_ID);
    });

    test(`Marker scores are always greater than or equal to ${FMG_GENE_STRENGTH_THRESHOLD}`, async ({
      page,
    }) => {
      await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

      clickDropdownOptionByName({
        page,
        testId: ADD_TISSUE_ID,
        name: "lung",
      });

      getCellTypeFmgButtonAndClick(page, "club cell");

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
    test("Display gene info panel in sidebar", async ({ page }) => {
      await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

      clickUntilOptionsShowUp({ page, testId: ADD_TISSUE_ID });
      selectFirstNOptions(1, page);

      clickUntilOptionsShowUp({ page, testId: ADD_GENE_ID });
      selectFirstNOptions(3, page);

      waitForHeatmapToRender(page);

      getFirstButtonAndClick(page, GENE_INFO_BUTTON_X_AXIS_TEST_ID);

      waitForElement(page, GENE_INFO_HEADER_TEST_ID);
      waitForElement(page, GENE_INFO_GENE_SYNONYMS_TEST_ID);
      waitForElement(page, GENE_INFO_NCBI_LINK_TEST_ID);
      waitForElement(page, RIGHT_SIDEBAR_TITLE_TEST_ID);

      getButtonAndClick(page, RIGHT_SIDEBAR_CLOSE_BUTTON_TEST_ID);

      waitForElementToBeRemoved(page, RIGHT_SIDEBAR_TITLE_TEST_ID);
    });

    test("Display gene info bottom drawer in cell info sidebar", async ({
      page,
    }) => {
      await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

      clickDropdownOptionByName({
        page,
        testId: ADD_TISSUE_ID,
        name: "lung",
      });

      clickUntilOptionsShowUp({ page, testId: ADD_GENE_ID });
      selectFirstNOptions(3, page);

      waitForHeatmapToRender(page);

      getCellTypeFmgButtonAndClick(page, "muscle cell");

      getFirstButtonAndClick(page, GENE_INFO_BUTTON_CELL_INFO_TEST_ID);

      waitForElement(page, GENE_INFO_HEADER_TEST_ID);
      waitForElement(page, GENE_INFO_GENE_SYNONYMS_TEST_ID);
      waitForElement(page, GENE_INFO_NCBI_LINK_TEST_ID);
      waitForElement(page, GENE_INFO_TITLE_SPLIT_TEST_ID);

      getButtonAndClick(page, GENE_INFO_CLOSE_BUTTON_SPLIT_TEST_ID);

      waitForElementToBeRemoved(page, GENE_INFO_TITLE_SPLIT_TEST_ID);

      getButtonAndClick(page, RIGHT_SIDEBAR_CLOSE_BUTTON_TEST_ID);

      waitForElementToBeRemoved(page, RIGHT_SIDEBAR_TITLE_TEST_ID);

      getFirstButtonAndClick(page, GENE_INFO_BUTTON_X_AXIS_TEST_ID);

      waitForElement(page, GENE_INFO_HEADER_TEST_ID);
      waitForElement(page, RIGHT_SIDEBAR_TITLE_TEST_ID);
    });
  });

  describe("Export CSV ", () => {
    test("Download CSV and validate length - no compare", async ({ page }) => {
      await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

      const SELECT_N_GENES = 3;

      // Select tissue
      clickUntilOptionsShowUp({ page, testId: ADD_TISSUE_ID });
      selectFirstOption(page);

      // Select gene
      clickUntilOptionsShowUp({ page, testId: ADD_GENE_ID });
      selectFirstNOptions(SELECT_N_GENES, page);

      waitForHeatmapToRender(page);

      clickUntilDownloadModalShowsUp({
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

      expect(fileName).toBe("blood.csv");

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
      expect(data.length).toBe(cellTypes.length * SELECT_N_GENES);
    });

    test("Download CSV and validate length - compare", async ({ page }) => {
      await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

      const SELECT_N_GENES = 3;

      // Select tissue
      clickUntilOptionsShowUp({ page, testId: ADD_TISSUE_ID });
      selectFirstOption(page);

      // Select gene
      clickUntilOptionsShowUp({ page, testId: ADD_GENE_ID });
      selectFirstNOptions(SELECT_N_GENES, page);

      // Select a compare dimension
      clickDropdownOptionByName({
        name: "Sex",
        page,
        testId: COMPARE_DROPDOWN_ID,
      });

      waitForHeatmapToRender(page);

      clickUntilDownloadModalShowsUp({
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

      expect(fileName).toBe("blood.csv");

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
      expect(data.length).toBe(cellTypes.length * SELECT_N_GENES);
    });
  });

  describe("Multiple file download - one tissue", () => {
    test("Download zip of all outputs (png,svg,csv) for one tissue", async ({
      page,
    }) => {
      await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

      clickUntilOptionsShowUp({ page, testId: ADD_TISSUE_ID });
      selectFirstOption(page);

      clickUntilOptionsShowUp({ page, testId: ADD_GENE_ID });
      selectFirstOption(page);

      waitForHeatmapToRender(page);

      clickUntilDownloadModalShowsUp({
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

      const files = ["blood.csv", "blood.png", "blood.svg"];

      expect(zipEntries.length).toBe(3);

      for (const entry of zipEntries) {
        expect(files.includes(entry.name)).toBe(true);
      }
    });
  });

  describe("Multiple file download - two tissues", () => {
    test("Download zip of all outputs (png,svg,csv) for two tissues", async ({
      page,
    }) => {
      await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

      // Select first tissue
      clickUntilOptionsShowUp({ page, testId: ADD_TISSUE_ID });
      selectFirstOption(page);

      // Select second tissue
      await clickUntilOptionsShowUp({ page, testId: ADD_TISSUE_ID });
      await selectNthOption(page, 2);

      clickUntilOptionsShowUp({ page, testId: ADD_GENE_ID });
      selectFirstOption(page);

      waitForHeatmapToRender(page);

      clickUntilDownloadModalShowsUp({
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
        "blood.csv",
        "blood.png",
        "blood.svg",
        "lung.csv",
        "lung.png",
        "lung.svg",
      ];

      expect(zipEntries.length).toBe(6);

      for (const entry of zipEntries) {
        expect(files.includes(entry.name)).toBe(true);
      }
    });
  });
});
