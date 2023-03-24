import { ElementHandle, expect, Page, test } from "@playwright/test";
import { ROUTES } from "src/common/constants/routes";
import type { RawPrimaryFilterDimensionsResponse } from "src/common/queries/wheresMyGene";
import { goToPage, isDevStagingProd, tryUntil } from "tests/utils/helpers";
import { TEST_URL } from "../common/constants";
import { getTestID, getText } from "../utils/selectors";
import { TISSUE_DENY_LIST } from "./fixtures/wheresMyGene/tissueRollup";
import fs from "fs";
import { parse } from "csv-parse/sync";
import AdmZip from "adm-zip";

const HOMO_SAPIENS_TERM_ID = "NCBITaxon:9606";

const GENE_LABELS_ID = "gene-labels";
const CELL_TYPE_LABELS_ID = "cell-type-name";
const ADD_TISSUE_ID = "add-tissue";
const ADD_GENE_ID = "add-gene";
const GENE_DELETE_BUTTON = "gene-delete-button";
const SOURCE_DATA_BUTTON_ID = "source-data-button";
const SOURCE_DATA_LIST_SELECTOR = `[data-test-id="source-data-list"]`;
const DOWNLOAD_BUTTON_ID = "download-button";

const MUI_CHIP_ROOT = ".MuiChip-root";
const FILTERS_PANEL_NOT_FOUND = "Filters panel not found";

const CELL_TYPE_SANITY_CHECK_NUMBER = 100;

const COMPARE_DROPDOWN_ID = "compare-dropdown";

const EXPORT_OUTPUT_DIR = "playwright-report/";

const { describe, skip } = test;

describe("Where's My Gene", () => {
  skip(!isDevStagingProd, "WMG BE API does not work locally or in rdev");

  test("renders the getting started UI", async ({ page }) => {
    await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

    // Getting Started section
    await expect(page).toHaveSelector(getText("STEP 1"));
    await expect(page).toHaveSelector(getText("Add Tissues"));

    await expect(page).toHaveSelector(getText("STEP 2"));
    await expect(page).toHaveSelector(getText("Add Genes"));

    await expect(page).toHaveSelector(getText("STEP 3"));
    await expect(page).toHaveSelector(getText("Explore Gene Expression"));

    // Filters Panel
    // (thuang): `*` is for intermediate match
    // https://playwright.dev/docs/selectors#intermediate-matches

    async function getTissueSelectorButton() {
      return page.$(getTestID(ADD_TISSUE_ID));
    }

    async function getGeneSelectorButton() {
      return page.$(getTestID(ADD_GENE_ID));
    }

    await clickUntilOptionsShowUp(getTissueSelectorButton, page);
    await selectFirstOption(page);

    await clickUntilOptionsShowUp(getGeneSelectorButton, page);
    await selectFirstOption(page);

    const filtersPanel = await page.$("*css=div >> text=Filters");

    await expect(filtersPanel).toHaveSelector(getText("Dataset"));
    await expect(filtersPanel).toHaveSelector(getText("Disease"));
    await expect(filtersPanel).toHaveSelector(
      getText("Self-Reported Ethnicity")
    );
    await expect(filtersPanel).toHaveSelector(getText("Sex"));

    // Legend
    const Legend = await page.$("*css=div >> text=Gene Expression");

    await expect(Legend).toHaveSelector(getText("Gene Expression"));
    await expect(Legend).toHaveSelector(getText("Expressed in Cells (%)"));

    // Info Panel
    const InfoPanel = await page.$("*css=div >> text=Filters");

    await expect(InfoPanel).toHaveSelector(getText("Methodology"));
    await expect(InfoPanel).toHaveSelector(
      getText("After filtering cells with low coverage ")
    );
    await expect(InfoPanel).toHaveSelector(getText("Source Data"));
  });

  test("Filters and Heatmap", async ({ page }) => {
    await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

    async function getTissueSelectorButton() {
      return page.$(getTestID(ADD_TISSUE_ID));
    }

    async function getGeneSelectorButton() {
      return page.$(getTestID(ADD_GENE_ID));
    }

    await clickUntilOptionsShowUp(getTissueSelectorButton, page);
    await selectFirstOption(page);

    await clickUntilOptionsShowUp(getGeneSelectorButton, page);
    await selectFirstOption(page);

    await waitForHeatmapToRender(page);

    const sexSelector = await getSexSelector();

    if (!sexSelector) throw Error("No sexSelector found");

    const selectedSexesBefore = await sexSelector.$$(MUI_CHIP_ROOT);

    await expect(selectedSexesBefore.length).toBe(0);

    await clickUntilOptionsShowUp(getSexSelectorButton, page);

    await selectFirstOption(page);

    const selectedSexesAfter = await sexSelector.$$(MUI_CHIP_ROOT);

    await expect(selectedSexesAfter.length).toBe(1);

    async function getFiltersPanel() {
      return page.$(getTestID("filters-panel"));
    }

    async function getSexSelector() {
      const filtersPanel = await getFiltersPanel();

      if (!filtersPanel) {
        throw Error(FILTERS_PANEL_NOT_FOUND);
      }

      return filtersPanel.$(getTestID("sex-filter"));
    }

    async function getSexSelectorButton() {
      const filtersPanel = await getFiltersPanel();

      if (!filtersPanel) {
        throw Error(FILTERS_PANEL_NOT_FOUND);
      }

      await filtersPanel.$(getTestID("sex-filter"));
      return filtersPanel.$("*css=button >> text=Sex");
    }
  });

  test("Primary and secondary filter crossfiltering", async ({ page }) => {
    await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

    const diseaseSelector = await getDiseaseSelector();
    const datasetSelector = await getDatasetSelector();

    if (!diseaseSelector) throw Error("No diseaseSelector found");
    if (!datasetSelector) throw Error("No datasetSelector found");

    await clickUntilOptionsShowUp(getTissueSelectorButton, page);
    const tissueOptionsBefore = await page.$$("[role=option]");
    const numberOfTissuesBefore = tissueOptionsBefore.length;
    await page.keyboard.press("Escape");

    await clickUntilOptionsShowUp(getDiseaseSelectorButton, page);
    const diseaseOptionsBefore = await page.$$("[role=option]");
    const numberOfDiseasesBefore = diseaseOptionsBefore.length;
    await page.keyboard.press("Escape");

    await clickUntilOptionsShowUp(getDatasetSelectorButton, page);
    const datasetOptions = await page.$$("[role=option]");
    await datasetOptions[0].click();
    await page.keyboard.press("Escape");

    await clickUntilOptionsShowUp(getTissueSelectorButton, page);
    const tissueOptionsAfter = await page.$$("[role=option]");
    const numberOfTissuesAfter = tissueOptionsAfter.length;
    await page.keyboard.press("Escape");

    await clickUntilOptionsShowUp(getDiseaseSelectorButton, page);
    const diseaseOptionsAfter = await page.$$("[role=option]");
    const numberOfDiseasesAfter = diseaseOptionsAfter.length;
    await page.keyboard.press("Escape");

    expect(numberOfDiseasesBefore).toBeGreaterThan(numberOfDiseasesAfter);
    expect(numberOfTissuesBefore).toBeGreaterThan(numberOfTissuesAfter);

    async function getTissueSelectorButton() {
      return page.$(getTestID(ADD_TISSUE_ID));
    }

    async function getFiltersPanel() {
      return page.$(getTestID("filters-panel"));
    }

    async function getDiseaseSelector() {
      const filtersPanel = await getFiltersPanel();

      if (!filtersPanel) {
        throw Error(FILTERS_PANEL_NOT_FOUND);
      }

      return filtersPanel.$("*css=div >> text=Disease");
    }

    async function getDiseaseSelectorButton() {
      const filtersPanel = await getFiltersPanel();

      if (!filtersPanel) {
        throw Error(FILTERS_PANEL_NOT_FOUND);
      }

      await filtersPanel.$("*css=div >> text=Disease");
      return filtersPanel.$("*css=button >> text=Disease");
    }

    async function getDatasetSelector() {
      const filtersPanel = await getFiltersPanel();

      if (!filtersPanel) {
        throw Error(FILTERS_PANEL_NOT_FOUND);
      }

      return filtersPanel.$("*css=div >> text=Dataset");
    }

    async function getDatasetSelectorButton() {
      const filtersPanel = await getFiltersPanel();

      if (!filtersPanel) {
        throw Error(FILTERS_PANEL_NOT_FOUND);
      }

      await filtersPanel.$("*css=div >> text=Dataset");
      return filtersPanel.$("*css=button >> text=Dataset");
    }
  });

  test("Selected tissue and no genes displays cell types", async ({ page }) => {
    await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

    async function getTissueSelectorButton() {
      return page.$(getTestID(ADD_TISSUE_ID));
    }

    await clickUntilOptionsShowUp(getTissueSelectorButton, page);
    await selectFirstOption(page);

    const cellTypes = await page.$$(getTestID("cell-type-name"));
    expect(cellTypes.length).toBeGreaterThan(0);
  });

  test("Source Data", async ({ page }) => {
    await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

    async function getSourceDataButton() {
      return page.$(getTestID(SOURCE_DATA_BUTTON_ID));
    }
    async function getTissueSelectorButton() {
      return page.$(getTestID(ADD_TISSUE_ID));
    }

    async function getGeneSelectorButton() {
      return page.$(getTestID(ADD_GENE_ID));
    }

    await clickUntilOptionsShowUp(getTissueSelectorButton, page);
    await selectFirstOption(page);

    await clickUntilOptionsShowUp(getGeneSelectorButton, page);
    await selectFirstOption(page);

    await waitForHeatmapToRender(page);
    await clickUntilSidebarShowsUp(getSourceDataButton, page);
    await expect(page).toHaveSelector(
      getText(
        "After filtering cells with low coverage (less than 500 genes expressed)"
      )
    );

    await tryUntil(
      async () => {
        const sourceDataList = await page.$(SOURCE_DATA_LIST_SELECTOR);

        if (!sourceDataList) throw Error("no source data displayed");

        const sourceDataListItems = await sourceDataList?.$$(
          ".MuiListItem-root"
        );

        expect(sourceDataListItems?.length).toBeGreaterThan(0);

        await page.mouse.click(0, 0);

        async function getFiltersPanel() {
          return page.$(getTestID("filters-panel"));
        }
        async function getDatasetSelector() {
          const filtersPanel = await getFiltersPanel();

          if (!filtersPanel) {
            throw Error(FILTERS_PANEL_NOT_FOUND);
          }

          return filtersPanel.$("*css=button >> text=Dataset");
        }
        const datasetSelector = await getDatasetSelector();
        if (!datasetSelector) throw Error("No datasetSelector found");
        const selectedDatasetsBefore = await datasetSelector.$$(MUI_CHIP_ROOT);
        await expect(selectedDatasetsBefore.length).toBe(0);
        await clickUntilOptionsShowUp(getDatasetSelector, page);
        await selectFirstOption(page);
        await clickUntilSidebarShowsUp(getSourceDataButton, page);

        const sourceDataListAfter = await page.$(SOURCE_DATA_LIST_SELECTOR);

        if (!sourceDataListAfter)
          throw Error(
            "no source data displayed after selecting dataset filter"
          );

        const sourceDataListAfterItems = await sourceDataListAfter?.$$(
          ".MuiListItem-root"
        );

        expect(sourceDataListAfterItems?.length).toBeGreaterThan(0);
      },
      { page }
    );
  });

  test("Hierarchical Clustering", async ({ page }) => {
    await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

    async function getTissueSelectorButton() {
      return page.$(getTestID(ADD_TISSUE_ID));
    }

    async function getGeneSelectorButton() {
      return page.$(getTestID(ADD_GENE_ID));
    }

    const TISSUE_COUNT = 1;
    const GENE_COUNT = 3;

    await clickUntilOptionsShowUp(getTissueSelectorButton, page);
    await selectFirstNOptions(TISSUE_COUNT, page);

    await clickUntilOptionsShowUp(getGeneSelectorButton, page);
    await selectFirstNOptions(GENE_COUNT, page);

    const beforeGeneNames = await getNames(
      `${getTestID(GENE_LABELS_ID)} span`,
      page
    );

    const beforeCellTypeNames = await getNames(
      `${getTestID(CELL_TYPE_LABELS_ID)}`,
      page
    );

    expect(beforeGeneNames.length).toBe(GENE_COUNT);

    // (thuang): Sometimes when API response is slow, we'll not capture all the
    // cell type names, so a sanity check that we expect at least 100 names
    expect(beforeCellTypeNames.length).toBeGreaterThan(
      CELL_TYPE_SANITY_CHECK_NUMBER
    );

    const cellTypeSortDropdown = await page.locator(
      getTestID("cell-type-sort-dropdown")
    );
    await cellTypeSortDropdown.click();
    await selectNthOption(2, page);

    const geneSortDropdown = await page.locator(
      getTestID("gene-sort-dropdown")
    );
    await geneSortDropdown.click();
    await selectNthOption(2, page);

    const afterGeneNames = await getNames(
      `${getTestID(GENE_LABELS_ID)} span`,
      page
    );

    const afterCellTypeNames = await getNames(
      `${getTestID(CELL_TYPE_LABELS_ID)}`,
      page
    );

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

  test("delete genes", async ({ page }) => {
    await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

    async function getTissueSelectorButton() {
      return page.$(getTestID(ADD_TISSUE_ID));
    }

    async function getGeneSelectorButton() {
      return page.$(getTestID(ADD_GENE_ID));
    }

    await clickUntilOptionsShowUp(getTissueSelectorButton, page);
    await selectFirstNOptions(1, page);

    await clickUntilOptionsShowUp(getGeneSelectorButton, page);
    await selectFirstNOptions(3, page);

    await waitForHeatmapToRender(page);

    const beforeGeneNames = await getNames(
      `${getTestID(GENE_LABELS_ID)} span`,
      page
    );

    await tryUntil(
      async () => {
        await page.keyboard.press("Backspace");

        // Testing single gene delete
        await page.hover(".gene-label-container");
        await page.click(getTestID(GENE_DELETE_BUTTON));

        const afterGeneNames = await getNames(
          `${getTestID(GENE_LABELS_ID)} span`,
          page
        );

        expect(afterGeneNames.length).toBe(beforeGeneNames.length - 1);
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

      async function getTissueSelectorButton() {
        return page.$(getTestID(ADD_TISSUE_ID));
      }

      async function getGeneSelectorButton() {
        return page.$(getTestID(ADD_GENE_ID));
      }

      await clickUntilOptionsShowUp(getTissueSelectorButton, page);
      await selectFirstNOptions(1, page);

      await clickUntilOptionsShowUp(getGeneSelectorButton, page);
      await selectFirstNOptions(3, page);

      await waitForHeatmapToRender(page);

      const beforeCellTypeNames = await getNames(
        `${getTestID(CELL_TYPE_LABELS_ID)}`,
        page
      );

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

      // Check all 3 Compare options work
      await clickDropdownOptionByName({
        name: "Disease",
        page,
        selector: getTestID(COMPARE_DROPDOWN_ID),
      });

      await tryUntil(
        async () => {
          const afterCellTypeNames = await getNames(
            `${getTestID(CELL_TYPE_LABELS_ID)}`,
            page
          );

          expect(
            afterCellTypeNames.find((name) => name.includes("  normal"))
          ).toBeTruthy();
        },
        { page }
      );

      await clickDropdownOptionByName({
        name: "Sex",
        page,
        selector: getTestID(COMPARE_DROPDOWN_ID),
      });

      await tryUntil(
        async () => {
          const afterCellTypeNames = await getNames(
            `${getTestID(CELL_TYPE_LABELS_ID)}`,
            page
          );

          expect(
            afterCellTypeNames.find((name) => name.includes("  female"))
          ).toBeTruthy();
        },
        { page }
      );

      await clickDropdownOptionByName({
        name: "Ethnicity",
        page,
        selector: getTestID(COMPARE_DROPDOWN_ID),
      });

      await tryUntil(
        async () => {
          const afterCellTypeNames = await getNames(
            `${getTestID(CELL_TYPE_LABELS_ID)}`,
            page
          );

          expect(
            afterCellTypeNames.find((name) => name.includes("  multiethnic"))
          ).toBeTruthy();
        },
        { page }
      );

      // Check selecting None removes stratification
      await clickDropdownOptionByName({
        name: "None",
        page,
        selector: getTestID(COMPARE_DROPDOWN_ID),
      });

      expect(
        beforeCellTypeNames.find((name) => name.includes("  "))
      ).toBeFalsy();
    });
  });

  describe("Export CSV ", () => {
    test("Download CSV and validate length - no compare", async ({ page }) => {
      await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

      const SELECT_N_GENES = 3;

      async function getTissueSelectorButton() {
        return page.$(getTestID(ADD_TISSUE_ID));
      }

      async function getGeneSelectorButton() {
        return page.$(getTestID(ADD_GENE_ID));
      }

      // Select tissue
      await clickUntilOptionsShowUp(getTissueSelectorButton, page);
      await selectFirstOption(page);

      // Select gene
      await clickUntilOptionsShowUp(getGeneSelectorButton, page);
      await selectFirstNOptions(SELECT_N_GENES, page);

      await waitForHeatmapToRender(page);

      const numCellTypes = (await page.$$(getTestID("cell-type-label-count")))
        .length;

      if (!numCellTypes) {
        throw Error("No cell types");
      }

      async function getDownloadButton() {
        return page.$(getTestID(DOWNLOAD_BUTTON_ID));
      }

      await clickUntilDownloadModalShowsUp(getDownloadButton, page);

      const pngCheckbox = await page.$(getTestID("png-checkbox"));
      if (!pngCheckbox) throw Error("No PNG checkbox");
      await pngCheckbox.click();

      const csvCheckbox = await page.$(getTestID("csv-checkbox"));
      if (!csvCheckbox) throw Error("No CSV checkbox");
      await csvCheckbox.click();

      // Start waiting for download before clicking. Note no await.
      const downloadPromise = page.waitForEvent("download");

      const downloadButton = await page.$(getTestID("dialog-download-button"));
      if (!downloadButton) throw Error("No download button");
      await downloadButton.click();

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
        from_line: 13, // Start on row with header names, skipping metadata
        skip_empty_lines: true,
      });

      // Validate number of rows
      expect(data.length).toBe(numCellTypes * SELECT_N_GENES);
    });

    test("Download CSV and validate length - compare", async ({ page }) => {
      await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

      const SELECT_N_GENES = 3;

      async function getTissueSelectorButton() {
        return page.$(getTestID(ADD_TISSUE_ID));
      }

      async function getGeneSelectorButton() {
        return page.$(getTestID(ADD_GENE_ID));
      }

      // Select tissue
      await clickUntilOptionsShowUp(getTissueSelectorButton, page);
      await selectFirstOption(page);

      // Select gene
      await clickUntilOptionsShowUp(getGeneSelectorButton, page);
      await selectFirstNOptions(SELECT_N_GENES, page);

      // Select a compare dimension
      await clickDropdownOptionByName({
        name: "Sex",
        page,
        selector: getTestID(COMPARE_DROPDOWN_ID),
      });

      await waitForHeatmapToRender(page);

      const numCellTypes = (await page.$$(getTestID("cell-type-label-count")))
        .length;

      if (!numCellTypes) {
        throw Error("No cell types");
      }

      async function getDownloadButton() {
        return page.$(getTestID(DOWNLOAD_BUTTON_ID));
      }

      await clickUntilDownloadModalShowsUp(getDownloadButton, page);

      const pngCheckbox = await page.$(getTestID("png-checkbox"));
      if (!pngCheckbox) throw Error("No PNG checkbox");
      await pngCheckbox.click();

      const csvCheckbox = await page.$(getTestID("csv-checkbox"));
      if (!csvCheckbox) throw Error("No CSV checkbox");
      await csvCheckbox.click();

      // Start waiting for download before clicking. Note no await.
      const downloadPromise = page.waitForEvent("download");

      const downloadButton = await page.$(getTestID("dialog-download-button"));
      if (!downloadButton) throw Error("No download button");
      await downloadButton.click();

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
        from_line: 13, // Start on row with header names, skipping metadata
        skip_empty_lines: true,
      });

      // Validate number of rows
      expect(data.length).toBe(numCellTypes * SELECT_N_GENES);
    });
  });

  describe("Multiple file download - one tissue", () => {
    test("Download zip of all outputs (png,svg,csv) for one tissue", async ({
      page,
    }) => {
      await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

      async function getTissueSelectorButton() {
        return page.$(getTestID(ADD_TISSUE_ID));
      }

      async function getGeneSelectorButton() {
        return page.$(getTestID(ADD_GENE_ID));
      }

      await clickUntilOptionsShowUp(getTissueSelectorButton, page);
      await selectFirstOption(page);

      await clickUntilOptionsShowUp(getGeneSelectorButton, page);
      await selectFirstOption(page);

      await waitForHeatmapToRender(page);

      async function getDownloadButton() {
        return page.$(getTestID(DOWNLOAD_BUTTON_ID));
      }

      await clickUntilDownloadModalShowsUp(getDownloadButton, page);

      const svgCheckbox = await page.$(getTestID("svg-checkbox"));
      if (!svgCheckbox) throw Error("No SVG checkbox");
      await svgCheckbox.click();

      const csvCheckbox = await page.$(getTestID("csv-checkbox"));
      if (!csvCheckbox) throw Error("No CSV checkbox");
      await csvCheckbox.click();

      // Start waiting for download before clicking. Note no await.
      const downloadPromise = page.waitForEvent("download");

      const downloadButton = await page.$(getTestID("dialog-download-button"));
      if (!downloadButton) throw Error("No download button");
      await downloadButton.click();

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

      async function getTissueSelectorButton() {
        return page.$(getTestID(ADD_TISSUE_ID));
      }

      async function getGeneSelectorButton() {
        return page.$(getTestID(ADD_GENE_ID));
      }

      // Select first tissue
      await clickUntilOptionsShowUp(getTissueSelectorButton, page);
      await selectFirstOption(page);

      // Select second tissue
      await clickUntilOptionsShowUp(getTissueSelectorButton, page);
      await selectNthOption(2, page);

      await clickUntilOptionsShowUp(getGeneSelectorButton, page);
      await selectFirstOption(page);

      await waitForHeatmapToRender(page);

      async function getDownloadButton() {
        return page.$(getTestID(DOWNLOAD_BUTTON_ID));
      }

      await clickUntilDownloadModalShowsUp(getDownloadButton, page);

      const svgCheckbox = await page.$(getTestID("svg-checkbox"));
      if (!svgCheckbox) throw Error("No SVG checkbox");
      await svgCheckbox.click();

      const csvCheckbox = await page.$(getTestID("csv-checkbox"));
      if (!csvCheckbox) throw Error("No CSV checkbox");
      await csvCheckbox.click();

      // Start waiting for download before clicking. Note no await.
      const downloadPromise = page.waitForEvent("download");

      const downloadButton = await page.$(getTestID("dialog-download-button"));
      if (!downloadButton) throw Error("No download button");
      await downloadButton.click();

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

async function getNames(selector: string, page: Page): Promise<string[]> {
  const geneLabelsLocator = await page.locator(selector);
  await tryUntil(
    async () => {
      const names = await geneLabelsLocator.allTextContents();
      expect(typeof names[0]).toBe("string");
    },
    { page }
  );

  return geneLabelsLocator.allTextContents();
}

async function clickUntilOptionsShowUp(
  getTarget: () => Promise<ElementHandle<SVGElement | HTMLElement> | null>,
  page: Page
) {
  await tryUntil(
    async () => {
      const target = await getTarget();

      if (!target) throw Error("no target");

      await target.click();
      const tooltip = await page.$("[role=tooltip]");

      if (!tooltip) throw Error("no tooltip");

      const options = await tooltip.$$("[role=option]");

      if (!options?.length) throw Error("no options");
    },
    { page }
  );
}

async function clickUntilDownloadModalShowsUp(
  getTarget: () => Promise<ElementHandle<SVGElement | HTMLElement> | null>,
  page: Page
) {
  await tryUntil(
    async () => {
      const target = await getTarget();

      if (!target) throw Error("no target");

      await target.click();

      const modal = await page.$("[class=bp4-dialog]");

      if (modal) throw Error("No download modal");
    },
    { page }
  );
}

async function clickUntilSidebarShowsUp(
  getTarget: () => Promise<ElementHandle<SVGElement | HTMLElement> | null>,
  page: Page
) {
  await tryUntil(
    async () => {
      const target = await getTarget();

      if (!target) throw Error("no target");

      await target.click();
      const drawer = await page.$("[class=bp4-drawer-header]");

      if (!drawer) throw Error("no drawer");
    },
    { page }
  );
}

// (thuang): This only works when a dropdown is open
async function selectFirstOption(page: Page) {
  await selectFirstNOptions(1, page);
}

async function selectFirstNOptions(count: number, page: Page) {
  for (let i = 0; i < count; i++) {
    await page.keyboard.press("ArrowDown");
    await page.keyboard.press("Enter");
  }

  await page.keyboard.press("Escape");
}

async function selectNthOption(number: number, page: Page) {
  // (thuang): Since the first option is now active, we need to offset by 1
  const step = number - 1;

  for (let i = 0; i < step; i++) {
    await page.keyboard.press("ArrowDown");
  }

  await page.keyboard.press("Enter");
  await page.keyboard.press("Escape");
}

async function waitForHeatmapToRender(page: Page) {
  await tryUntil(
    async () => {
      const canvases = await page.$$("canvas");
      await expect(canvases.length).not.toBe(0);
    },
    { page }
  );
}

async function clickDropdownOptionByName({
  page,
  selector,
  name,
}: {
  page: Page;
  selector: string;
  name: string;
}) {
  const dropdown = await page.locator(selector);
  await dropdown.click();

  const option = await page.locator(`[role=option] >> text=${name}`);
  await option.click();
}
