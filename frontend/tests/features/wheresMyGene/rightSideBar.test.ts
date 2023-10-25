/**
 * Test suite for select filter-related utils.
 */
import { expect } from "@playwright/test";
import { goToWMG } from "../../utils/wmgUtils";
import { selectFirstOption } from "tests/utils/helpers";
import { getById } from "tests/utils/selectors";
import { test } from "tests/common/test";

const { describe } = test;

describe("Right side bar", () => {
  test("Should link out to cellxgene documentation", async ({
    page,
    context,
  }) => {
    //click  documentation under Help and documentation this opens on a new page
    // navigate to url
    await goToWMG(page);
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

  test("Should load right side bar loads with legend, methodology, and source data", async ({
    page,
  }) => {
    // navigate to url
    await goToWMG(page);
    //verify gene-expression color scale is visible
    await expect(
      page.locator(getById("visualization-color-scale"))
    ).toBeVisible();

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
  test("Should links out to scExpression documentation", async ({
    page,
    context,
  }) => {
    // navigate to url
    await goToWMG(page);
    //click the source data icon
    await page.getByTestId("source-data-button").click();
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
});
