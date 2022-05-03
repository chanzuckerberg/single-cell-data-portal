import { ElementHandle } from "playwright";
import { ROUTES } from "src/common/constants/routes";
import { goToPage, tryUntil } from "tests/utils/helpers";
import { TEST_ENV, TEST_URL } from "../common/constants";
import { getTestID, getText } from "../utils/selectors";

//(thuang): BE API doesn't work in local happy
const TEST_ENVS = ["dev", "staging", "prod"];

const describeIfDevStagingProd = TEST_ENVS.includes(TEST_ENV)
  ? describe
  : describe.skip;

describeIfDevStagingProd("Where's My Gene", () => {
  it("renders the getting started UI", async () => {
    await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`);

    // Getting Started section
    await expect(page).toHaveSelector(getText("Getting Started"));
    await expect(page).toHaveSelector(
      getText(
        "Visualize and compare the expression of genes across cell types using a dot plot"
      )
    );

    // Step 1
    await expect(page).toHaveSelector(
      getText("Add tissues you are interested in exploring")
    );

    // Step 2
    await expect(page).toHaveSelector(
      getText("Add genes of interest to the visualization")
    );

    // Step 3
    await expect(page).toHaveSelector(
      getText("Darker dots represent higher relative gene expression")
    );

    // Beta callout
    await expect(page).toHaveSelector(getText("This feature is in beta"));

    // Filters Panel
    // (thuang): `*` is for intermediate match
    // https://playwright.dev/docs/selectors#intermediate-matches
    const filtersPanel = await page.$("*css=div >> text=Filters");

    await expect(filtersPanel).toHaveSelector(getText("Dataset"));
    await expect(filtersPanel).toHaveSelector(getText("Disease"));
    await expect(filtersPanel).toHaveSelector(getText("Ethnicity"));
    await expect(filtersPanel).toHaveSelector(getText("Sex"));

    // Info Panel
    const InfoPanel = await page.$("*css=div >> text=Info");

    await expect(InfoPanel).toHaveSelector(getText("Gene Expression"));
    await expect(InfoPanel).toHaveSelector(getText("Expressed in Cells (%)"));
    await expect(InfoPanel).toHaveSelector(getText("Methodology"));
    await expect(InfoPanel).toHaveSelector(
      getText("After filtering cells with low coverage ")
    );
    await expect(InfoPanel).toHaveSelector(getText("Source Data"));
  });

  test("Filters and Heatmap", async () => {
    await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`);

    async function getTissueSelectorButton() {
      return page.$(getTestID("add-tissue"));
    }

    async function getGeneSelectorButton() {
      return page.$(getTestID("add-gene"));
    }

    await clickUntilOptionsShowUp(getTissueSelectorButton);
    await selectFirstOption();

    await clickUntilOptionsShowUp(getGeneSelectorButton);
    await selectFirstOption();

    await tryUntil(async () => {
      const canvases = await page.$$("canvas");
      await expect(canvases.length).not.toBe(0);
    });

    const sexSelector = await getSexSelector();

    if (!sexSelector) throw Error("No sexSelector found");

    const selectedSexesBefore = await sexSelector.$$(".MuiChip-root");

    await expect(selectedSexesBefore.length).toBe(0);

    await clickUntilOptionsShowUp(getSexSelectorButton);

    await selectFirstOption();

    const selectedSexesAfter = await sexSelector.$$(".MuiChip-root");

    await expect(selectedSexesAfter.length).toBe(1);

    async function getFiltersPanel() {
      return page.$(getTestID("filters-panel"));
    }

    async function getSexSelector() {
      const filtersPanel = await getFiltersPanel();

      if (!filtersPanel) {
        throw Error("Filters panel not found");
      }

      return filtersPanel.$("*css=div >> text=Sex");
    }

    async function getSexSelectorButton() {
      const filtersPanel = await getFiltersPanel();

      if (!filtersPanel) {
        throw Error("Filters panel not found");
      }

      await filtersPanel.$("*css=div >> text=Sex");
      return filtersPanel.$("*css=button >> text=Sex");
    }
  });
});

async function clickUntilOptionsShowUp(
  getTarget: () => Promise<ElementHandle<SVGElement | HTMLElement> | null>
) {
  await tryUntil(async () => {
    const target = await getTarget();

    if (!target) throw Error("no target");

    await target.click();
    const tooltip = await page.$("[role=tooltip]");

    if (!tooltip) throw Error("no tooltip");

    const options = await tooltip.$$("[role=option]");

    if (!options?.length) throw Error("no options");
  });
}

// (thuang): This only works when a dropdown is open
async function selectFirstOption() {
  await page.keyboard.press("ArrowDown");
  await page.keyboard.press("Enter");
  await page.keyboard.press("Escape");
}
