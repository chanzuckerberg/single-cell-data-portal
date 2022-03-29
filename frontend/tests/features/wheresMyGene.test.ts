import { ElementHandle } from "playwright";
import { ROUTES } from "src/common/constants/routes";
import { goToPage, tryUntil } from "tests/utils/helpers";
import { TEST_URL } from "../common/constants";
import { getTestID, getText } from "../utils/selectors";

// DEBUG
// DEBUG
// DEBUG
// Temporarily skip WMG tests until BE API is available on deployed envs
// const TEST_ENVS = ["dev", "staging"];

// const describeIfDevOrStaging = TEST_ENVS.includes(TEST_ENV)
//   ? describe
//   : describe.skip;

// describeIfDevOrStaging("Where's My Gene", () => {
describe.skip("Where's My Gene", () => {
  it("renders the expected elements", async () => {
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
    await expect(filtersPanel).toHaveSelector(getText("Development Stage"));
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

  test("Filters", async () => {
    await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`);

    const filtersPanel = await page.$("*css=div >> text=Filters");

    if (!filtersPanel) {
      throw Error("Filters panel not found");
    }

    const sexSelector = await filtersPanel.$("*css=div >> text=Sex");

    const sexSelectorButton = await filtersPanel.$("*css=button >> text=Sex");

    if (!sexSelector || !sexSelectorButton) {
      throw Error("Dataset selector or button not found");
    }

    const selectedSexesBefore = await sexSelector.$$(".MuiChip-root");

    await expect(selectedSexesBefore.length).toBe(0);

    await clickUntilOptionsShowUp(sexSelectorButton);

    await selectFirstOption();

    const selectedSexesAfter = await sexSelector.$$(".MuiChip-root");

    await expect(selectedSexesAfter.length).toBe(1);
  });

  test("Heatmap", async () => {
    await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`);

    const geneSelectorButton = await page.$(getTestID("add-gene"));

    if (!geneSelectorButton) {
      throw Error("Gene selector button not found");
    }

    await clickUntilOptionsShowUp(geneSelectorButton);

    await selectFirstOption();

    await tryUntil(async () => {
      const canvases = await page.$$("canvas");
      await expect(canvases.length).not.toBe(0);
    });
  });
});

async function clickUntilOptionsShowUp(
  target: ElementHandle<SVGElement | HTMLElement> | null
) {
  await tryUntil(async () => {
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
