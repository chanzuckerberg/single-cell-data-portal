/**
 * Test suite for select filter-related utils.
 */
import { expect, test } from "@playwright/test";
import { goToWMG } from "../utils/wmgUtils";
import { isDevStagingProd } from "tests/utils/helpers";
import { getID, getTestID, getText } from "tests/utils/selectors";
import { selectFirstOption } from "./wheresMyGene.test";

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
      await page.click(getTestID("InputDropdown")),
      await selectFirstOption(page),
    ]);

    // wait until the new page fully loads
    await newPage.waitForLoadState();

    // expect the header on the new page to be visible
    expect(
      newPage.locator(
        getID("gene-expression--query-gene-expression-across-tissues")
      )
    ).toBeVisible();
  });

  test("Clicking on cellxgene logo takes you back to data portal", async ({
    page,
  }) => {
    //click the  logo on the top bar opens  opens on a new page
    await Promise.all([
      page.waitForNavigation({ waitUntil: "load" }),
      page.click(getTestID("logo")),
    ]);

    // expect the header on the new page to be visible
    expect(page.locator("h1")).toContainText(
      "Discover the mechanisms of human health"
    );
  });
  test.only("Right side bar loads with legend, methodology, and source data", async ({
    page,
  }) => {
    //verify gene-expression color scale is visible
    await expect(
      page.locator(getID("visualization-color-scale"))
    ).toBeVisible();
    await page.locator(getText("Loading")).waitFor({ state: "visible" });

    //click the source data icon
    await page.getByTestId("source-data-button").click();

    // expect the header to be visible
    await expect(page.locator("h4")).toContainText("Source Data");
    // expect source data to load
    await expect(
      page.locator('[data-testid="source-data-list"] a').count()
    ).toBeGreaterThan(0);
    // expect the header to be visible
    await expect(page.locator("h5")).toContainText("Methodology");
  });
  test("Link in methodology section links out to scExpression documentation", async ({
    page,
    context,
  }) => {
    //click the source data icon
    await page.getByTestId('[data-testid="source-data-button"]').click();

    const [newPage] = await Promise.all([
      context.waitForEvent("page"),
      //click our documentation link
      page.locator("a").locator(getText("our documentation")),
    ]);

    // wait until the new page fully loads
    await newPage.waitForLoadState();

    // expect the header on the new page to be visible
    expect(
      newPage.locator(
        getID("gene-expression--query-gene-expression-across-tissues")
      )
    ).toBeVisible();
  });
});
