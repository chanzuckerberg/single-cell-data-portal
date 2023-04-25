/**
 * Test suite for select filter-related utils.
 */
import { expect, test } from "@playwright/test";
import { goToWMG } from "../../utils/wmgUtils";
import { isDevStagingProd, selectFirstOption } from "tests/utils/helpers";
import { getById } from "tests/utils/selectors";

const { describe, skip } = test;

describe("Right side bar", () => {
  skip(!isDevStagingProd, "WMG BE API does not work locally or in rdev");
  test.beforeEach(async ({ page }) => {
    // navigate to url
    await goToWMG(page);
  });

  test("Documentation should link out to cellxgene documentation", async ({
    page,
    context,
  }) => {
    //click  documentation under Help and documentation this opens on a new page
    const [newPage] = await Promise.all([
      context.waitForEvent("page"),
      await page.getByTestId("InputDropdown").click(),
      await selectFirstOption(page),
    ]);

    // wait until the new page to fully loads
    await newPage.waitForLoadState();

    // expect the header on the new page to be visible
    expect(
      newPage.locator(
        getById("gene-expression--query-gene-expression-across-tissues")
      )
    ).toBeVisible();
  });

  test("Clicking on cellxgene logo takes you back to data portal", async ({
    page,
  }) => {
    //click the  logo on the top bar opens  opens on a new page
    await Promise.all([
      page.waitForNavigation({ waitUntil: "load" }),
      page.getByTestId("logo").click(),
    ]);

    // expect the header on the new page to be visible
    expect(page.locator("h1")).toContainText(
      "Discover the mechanisms of human health"
    );
  });
  test("Right side bar loads with legend, methodology, and source data", async ({
    page,
  }) => {
    //verify gene-expression color scale is visible
    await expect(
      page.locator(getById("visualization-color-scale"))
    ).toBeVisible();
    await page.getByTestId("add-tissue-btn").click();

    //click the source data icon
    await page.getByTestId("source-data-button").click();

    // expect the header to be visible
    await expect(page.locator("h4")).toContainText("Source Data");

    // expect source data to load
    expect(
      await page.getByTestId("source-data-list").locator("a").count()
    ).toBeGreaterThan(0);

    // expect the header to be visible
    await expect(page.locator("h5")).toContainText("Methodology");
  });
  test("Link in methodology section links out to scExpression documentation", async ({
    page,
    context,
  }) => {
    //click the source data icon
    await page.getByTestId("source-data-button").click();

    //close notification popup
    await page.getByTestId("notificationCloseButton").click();

    const [newPage] = await Promise.all([
      context.waitForEvent("page"),
      //click our documentation link
      await page.locator("a").getByText("our documentation").click(),
    ]);

    // wait until the new page fully loads
    await newPage.waitForLoadState();

    // expect the header on the new page to be visible
    expect(
      newPage.locator(
        getById("gene-expression--query-gene-expression-across-tissues")
      )
    ).toBeVisible();
  });

  test("When scaled check box is checked, gene expression should scale from 0 to 1", async ({
    page,
  }) => {
    const COLOR_SCALE = '[id="relative-gene-expression"] .low-high';

    //verify the correct values are displayed on the color scale
    await expect(page.locator(COLOR_SCALE)).toContainText("0.0");
    await expect(page.locator(COLOR_SCALE)).toContainText("1.0");

    //click the color scale drop-down
    await page.getByTestId("color-scale-dropdown").click();

    //select the second option
    await page.locator('[data-option-index="1"]').click();

    //verify the correct values are displayed on the color scale
    await expect(page.locator(COLOR_SCALE)).toContainText("0.0");
    await expect(page.locator(COLOR_SCALE)).toContainText("6.0");
  });
});
